'''
Created on Sep 12, 2011

@author: uemit.seren
'''
import gviz_api, os
import h5py
import tables
import numpy,math
from datetime import datetime
from JBrowseDataSource import DataSource as JBrowseDataSource
import cherrypy
import cPickle
from rnaseq_records import RNASeqRecords
from gviz_api import *
from cherrypy.lib.static import serve_file

@staticmethod
def SingleValueToJSPerf(value, value_type, escape_func=None):
    if escape_func is None:
        escape_func = gviz_api.DataTable._EscapeValue
    if isinstance(value, tuple):
        # In case of a tuple, we run the same function on the value itself and
        # add the formatted value.
        if (len(value) not in [2, 3] or
          (len(value) == 3 and not isinstance(value[2], dict))):
            raise DataTableException("Wrong format for value and formatting - %s." %
                                 str(value))
        if not isinstance(value[1], types.StringTypes + (types.NoneType,)):
            raise DataTableException("Formatted value is not string, given %s." %
                                 type(value[1]))
        js_value = DataTable.SingleValueToJS(value[0], value_type)
        if value[1] is None:
            return (js_value, None)
        return (js_value, escape_func(value[1]))

    # The standard case - no formatting.
    t_value = type(value)
    if value is None:
        return "null"
    if value_type == "boolean":
        if value:
            return "true"
        return "false"
    elif value_type == "number":
        if isinstance(value, (int, long, float)):
            if isinstance(value,float):
                return "%.2f" % value
            else:
                return str(value)
        raise DataTableException("Wrong type %s when expected number" % t_value)

    elif value_type == "string":
        if isinstance(value, tuple):
            raise DataTableException("Tuple is not allowed as string value.")
        return escape_func(value)

    elif value_type == "date":
        if not isinstance(value, (datetime.date, datetime.datetime)):
            raise DataTableException("Wrong type %s when expected date" % t_value)
        # We need to shift the month by 1 to match JS Date format
        return "new Date(%d,%d,%d)" % (value.year, value.month - 1, value.day)

    elif value_type == "timeofday":
        if not isinstance(value, (datetime.time, datetime.datetime)):
            raise DataTableException("Wrong type %s when expected time" % t_value)
        return "[%d,%d,%d]" % (value.hour, value.minute, value.second)

    elif value_type == "datetime":
        if not isinstance(value, datetime.datetime):
            raise DataTableException("Wrong type %s when expected datetime" %
                                 t_value)
        return "new Date(%d,%d,%d,%d,%d,%d)" % (value.year,
                                              value.month - 1,  # To match JS
                                              value.day,
                                              value.hour,
                                              value.minute,
                                              value.second)
    # If we got here, it means the given value_type was not one of the
    # supported types.
    raise DataTableException("Unsupported type %s" % value_type)

gviz_api.DataTable.SingleValueToJS =  SingleValueToJSPerf

class RNASeqService:
    base_path = "/srv/data/rnaseq-web/"
    base_path_jbrowse = base_path
    base_path_datasets = base_path + "datasets/"
    base_jbrowse_path = base_path_jbrowse + "jbrowse_1.2.1/"
    track_folder = "TAIR10"

    #tracks/Chr%s/TAIR10/"
    __datasource = None
    _lazyArrayChunks = [{}, {}, {}, {}, {}]
    hdf5_filename = base_path + "rnaseq_data.hdf5"
    genomeStats_hdf5_filename = base_path+ 'genestats.hdf5'
    gene_annot_file = base_path + "genome_annotation.pickled"
    phenotype_file = base_path + "phenotypes_fake.pickled"

    def __init__(self):
        self.rnaseq_records = RNASeqRecords(self.hdf5_filename)
        self.accessions,self.accessions_ix = self.rnaseq_records.getAccessions()
        #self.accessions_ix = dict((p['ecotype'], i) for i, p in enumerate(people))
        gene_annot_file = open(self.gene_annot_file, 'rb')
        self.gene_annotation = cPickle.load(gene_annot_file)
        gene_annot_file.close()
        phenotype_file = open(self.phenotype_file,'rb')
        #self.phenotypes = cPickle.load(phenotype_file)
        #phenotype_file.close()
        self.genomestats_file = h5py.File(self.genomeStats_hdf5_filename,'r')
        self.genome_wide_stats =   [{'name':'genecount','label':'# Genes','isStackable':False,'isStepPlot':True}, \
                    {'name':'Dn'},{'name':'Ds'},{'name':'Pn'},{'name':'Ps'},{'name':'MK_test'}, \
                    {'name':'Pn/Ps'},{'name':'Dn/Ds'},{'name':'alpha'},{'name':'pesudo_tajima_d_s','label':'Tajima D (synonymous)'},{'name':'pesudo_tajima_d_ns','label':'Tajima D (non-synonymous)'}]


    def _getLocationDistribution(self):
        import bisect, itertools
        from operator import itemgetter
        locations = {}
        selector = lambda row:row['country']
        for country, rows in itertools.groupby(sorted(self.accessions, key=itemgetter('country')), selector):
            locations[country] = len(list(rows))
        return locations

    def getLocationDistributionData(self):
        location_dist = self._getLocationDistribution()
        data = []
        for country, count in location_dist.iteritems():
            data.append({'Country':country, 'Count':count})
        column_name_type_ls = [("Country", ("string", "Country")), \
                          ("Count", ("number", "Count"))]
        column_ls = [row[0] for row in column_name_type_ls]
        data_table = gviz_api.DataTable(dict(column_name_type_ls))
        data_table.LoadData(data)
        return data_table.ToJSon(columns_order=column_ls)

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getTopResults(self,environ_type,result_type,range_start=0,range_length=1000, gene='',snp_chr='',snp_pos ='',gene_chr = '',min_score=''):
        results, count,range_start = self.rnaseq_records.getTopResults(environ_type,result_type,int(range_start), int(range_length), gene, snp_chr,snp_pos,gene_chr,min_score)
        return {'results':results, 'count':count, 'start':int(range_start),'length':int(range_length)}



    @cherrypy.expose
    @cherrypy.tools.response_headers(headers=[('Content-Type','application/csv')])
    def downloadAssociationData(self,phenotype,environment,dataset=None,transformation=None,result=None):
        import StringIO,csv
        results = self.rnaseq_records.get_results_for_csv(phenotype,environment,dataset,transformation,result)
        result_name = result
        if environment != 'GxE' and result not in ['LM','KW','EX']:
            result_name = 'EX_step%s' % result[-1]
        filename = '%s_%s_%s.csv' % (phenotype,environment,result_name)
        cherrypy.response.headers['Content-disposition'] = 'attachment;filename=%s' %filename
        content = ''
        tempfile = StringIO.StringIO()
        writer = csv.writer(tempfile,delimiter=',')
        writer.writerows(results)
        content = tempfile.getvalue()
        tempfile.close()
        return content;


    """"""
    @cherrypy.expose
    def downloadTopResults(self,environ_type,result_type, gene='',snp_chr='',snp_pos ='',gene_chr = '',min_score=''):
        filename ="%s_%s" % (environ_type,result_type)
        path = self.base_path + "/%s.tar.gz" % filename
        if os.path.isfile(path):
            return serve_file(path,"application/gzip","attachment")
        else:
            path = self.base_path+"/%s.csv" %filename
        return serve_file(path, "application/csv", "attachment")

    def _getPhenotypes(self, range_start,range_length,name='', chr='', start='', end=''):
        return self.rnaseq_records.getPhenotypes(range_start, range_length, name, chr, start, end)

        """ OLD IMPLEMENTATION
        phenotypes_to_filter = []
        phenotypes_to_return = []
        stop = int(range_start)+int(range_length)

        isCriteria = False
        if (name != '') or (chr !='') or (start != '') or (end != ''):
            isCriteria = True
        if isCriteria == False:
            phenotypes_to_filter = self.phenotypes[:]
        else:
            i = 0
            for i in range(0,len(self.phenotypes)):
                phenotype = self.phenotypes[i]
                if (name =='' or name.lower() in unicode(phenotype['name'],'utf8').lower()) and \
                    (chr =='' or int(chr) == phenotype['chr']) and (start == '' or int(start) <= phenotype['start']) and \
                    (end == '' or int(end) >= phenotype['end']):
                    phenotypes_to_filter.append(phenotype)

        count = len(phenotypes_to_filter)
        if stop > count:
            range_start = range_start - stop
            stop = count
        if range_start < 0:
            range_start = 0
        for i in range(range_start,stop):
            phenotypes_to_return.append(phenotypes_to_filter[i])
        return phenotypes_to_return,count,range_start
        """

    def _getHistogramDataTable(self,phenotype,environment,dataset='Fullset',transformation='raw'):
        column_name_type_ls = [("x_axis", ("string","Phenotype Value")), \
                            ("frequency",("number", "Frequency"))]
        description = dict(column_name_type_ls)
        data_table = gviz_api.DataTable(description)
        data_table.LoadData(self.rnaseq_records.get_phenotype_bins(phenotype,environment,dataset,transformation))
        column_ls = [row[0] for row in column_name_type_ls]
        return data_table.ToJSon(columns_order=column_ls)

    def _getCombinedHistogramDataTable(self,phenotype,dataset='Fullset',transformation='raw'):
        column_name_type_ls = [("x_axis", ("string","Phenotype Value")), \
                            ("frequency",("number", "10 C ")),
                            ("frequency_alt",("number","16 C"))]
        description = dict(column_name_type_ls)
        data_table = gviz_api.DataTable(description)
        data_table.LoadData(self.rnaseq_records.get_phenotype_bins(phenotype,'both',dataset,transformation))
        column_ls = [row[0] for row in column_name_type_ls]
        return data_table.ToJSon(columns_order=column_ls)

    def _getPhenotypeExplorerData(self,phenotype,environment,dataset,transformation,phen_vals = None):
        import datetime
        result = {}
        column_name_type_ls = [("label", ("string","ID Name Phenotype")),("date", ("date", "Date")),\
                               ("lon",("number", "Longitude")),  ("lat",("number", "Latitude")), \
                               ("phenotype", ("number", "Phenotype")),("name", ("string", "Native Name")),\
                               ("country", ("string", "Country")) ]
        if phen_vals == None:
            phen_vals = self.rnaseq_records.get_phenotype_values(phenotype, environment, dataset, transformation)
        data = []
        for i, ecotype in enumerate(phen_vals['ecotype']):
            accession = self.accessions[self.accessions_ix.get(ecotype,-1)]
            if accession is not None:
                value = phen_vals['value'][i]
                label = '%s ID:%s Phenotype:%s.'%(accession['name'], ecotype, value)
                data.append({'label':label,'date':datetime.date(2009,2,3),'accession_id':ecotype,\
                             'lon':accession['longitude'],'lat':accession['latitude'],\
                             'phenotype':value,'name':accession['name'],'country':accession['country']})
        data_table = gviz_api.DataTable(dict(column_name_type_ls))
        data_table.LoadData(data)
        column_ls = [row[0] for row in column_name_type_ls]
        result= data_table.ToJSon(columns_order=column_ls)
        return result

    def _getGeneCountHistogramData(self,chr):
        if self.__datasource == None:
            self.__datasource = JBrowseDataSource(self.base_jbrowse_path,self.track_folder)
        bpPerBin,histogramData = self.__datasource.getGeneHistogram(chr)
        binCount = len(histogramData)
        maxValue = max(histogramData)
        positions = [i*bpPerBin for i in range(0,binCount+1)]
        return zip(positions,histogramData)

    @cherrypy.expose
    def getAccessions(self):
        return {'accessions':self.accessions}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getPhenotypes(self, range_start=0, range_length=1000, name='', chr='', start='', end=''):
        phenotypes, count,range_start = self._getPhenotypes(int(range_start), int(range_length), name, chr, start, end)
        return {'phenotypes':phenotypes, 'count':count, 'start':int(range_start),'length':int(range_length)}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getPhenotypeData(self,id):
        retval = {}
        phenotypes = self._getPhenotypes(0,1,id)
        phenotype_info = phenotypes[0][0]
        phenotype = self.rnaseq_records.get_phenotype_info(id)
        phenotype['start'] = phenotype_info['start']
        phenotype['end'] = phenotype_info['end']
        phenotype['chr'] = phenotype_info['chr']
        phenotype['phenotype_id'] = phenotype_info['phenotype_id']
        phenotype['maxScore10C'] = phenotype_info['maxScore10C']
        phenotype['maxScore16C'] = phenotype_info['maxScore16C']
        phenotype['maxScoreFull'] = phenotype_info['maxScoreFull']

        phenotype['pseudoHeritability10C'] = phenotype_info['pseudoHeritability10C']
        phenotype['pseudoHeritability16C'] = phenotype_info['pseudoHeritability16C']

        retval['phenotype'] = phenotype;
        retval['histogramdataTable'] = self._getCombinedHistogramDataTable(id)
        return retval

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getPhenotypeMotionChartData(self,phenotype,environment,dataset,transformation):
        return {"data":self._getPhenotypeExplorerData(phenotype,environment,dataset,transformation)}


    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getGxEGWASData(self,phenotype,result):
        try:
            import math
            retval = {}
            result =  result.lower()
            association_result = self.rnaseq_records.get_gxe_results_by_chromosome(phenotype,result)
            if result != 'combined':
                description = [('position',"number", "Position"),('value','number','-log Pvalue')]
                columns_order=("position", "value")
            else:
                description =  [('position',"number", "Position"),('full_value','number','Full model'),('genetic_value','number','Common effect'),('environ_value','number','Interaction effect')]
                columns_order=("position", "full_value","genetic_value","environ_value")
            chr2data ={}
            for i in range(1,6):
                if result == 'combined':
                    data = zip(association_result[i]['position'],association_result[i]['full'],association_result[i]['genetic'],association_result[i]['environ'])
                else:
                    data = zip(association_result[i]['position'].tolist(),association_result[i]['score'].tolist())
                data_table = gviz_api.DataTable(description)
                data_table.LoadData(data)
                chr2data[i] =  data_table.ToJSon()
            retval['chr2data'] = chr2data
            retval['chr2length'] = association_result['chromosome_ends']
            retval['max_value'] = association_result['max_score']
            if 'no_of_tests' in association_result:
                retval['bonferroniThreshold'] = -math.log10(1.0/(association_result['no_of_tests']*20.0))
            else:
                retval['bonferroniThreshold'] = -math.log10(1.0/(16000000.0*20.0))
        except Exception, err:
            retval ={"status":"ERROR","statustext":"%s"%str(err)}
        return retval


    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getGWASData(self,phenotype,environment,dataset,transformation,result):
        try:
            import math
            retval = {}
            association_result = self.rnaseq_records.get_results_by_chromosome(phenotype,environment,dataset,transformation,result)

            description = [('position',"number", "Position"),('value','number','-log Pvalue')]
            chr2data ={}
            for i in range(1,6):
                data = zip(association_result[i]['position'].tolist(),association_result[i]['score'].tolist())
                data.sort()
                data_table = gviz_api.DataTable(description)
                data_table.LoadData(data)
                chr2data[i] =  data_table.ToJSon(columns_order=("position", "value"))
            retval['chr2data'] = chr2data
            retval['chr2length'] = association_result['chromosome_ends']
            retval['max_value'] = association_result['max_score']
            if 'no_of_tests' in association_result:
                retval['bonferroniThreshold'] = -math.log10(1.0/(association_result['no_of_tests']*20.0))
            else:
                retval['bonferroniThreshold'] = -math.log10(1.0/(55000000.0*20.0))
        except Exception, err:
            retval ={"status":"ERROR","statustext":"%s"%str(err)}
        return retval

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getGenes(self,chromosome,start,end,isFeatures=''):
        try:
            if self.__datasource == None:
                self.__datasource = JBrowseDataSource(self.base_jbrowse_path,self.track_folder)
            genes = []
            genes = self.__datasource.getGenes(chromosome, int(start), int(end), bool(isFeatures))
            retval = {'status': 'OK','genes':genes}
        except Exception,err:
            retval =  {"status":"ERROR","statustext":"%s" %str(err)}
        return retval

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getGeneFromName(self,query):
        try:
            if self.__datasource == None:
                self.__datasource = JBrowseDataSource(self.base_jbrowse_path,self.track_folder)
            genes = []
            gene = self.__datasource.getGeneFromName(query)
            retval = {'status': 'OK','gene':gene}
        except Exception,err:
            retval =  {"status":"ERROR","statustext":"%s" %str(err)}
        return retval

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getGenesFromQuery(self,query):
        try:
            if self.__datasource == None:
                self.__datasource = JBrowseDataSource(self.base_jbrowse_path,self.track_folder)
            genes = []
            genes = self.__datasource.getGenesFromQuery(query)
            isMore = False
            count=0
            if len(genes) > 20:
                count = len(genes)
                genes = genes[0:20]
                isMore = True
            retval = {'status': 'OK','isMore':isMore,'count':count,'genes':genes}
        except Exception,err:
            retval =  {"status":"ERROR","statustext":"%s" %str(err)}
        return retval

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getGeneDescription(self,gene):
        try:
            gene_parts = gene.split('.')
            description = self.gene_annotation[gene_parts[0]][gene]['functional_description']['computational_description']
            retval = {'status': 'OK','description':description}
        except Exception,err:
            retval =  {"status":"ERROR","statustext":"%s" %str(err)}
        return retval

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getGenomeStatsList(self):
        try:
            retval = {'status': 'OK','stats':self.genome_wide_stats}
        except Exception,err:
            retval =  {"status":"ERROR","statustext":"%s" %str(err)}
        return retval

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getGenomeStatsData(self,stats,chr):
        try:
            if chr == "null":
                chr = "Chr1"
            description = [('position',"number", "Position")]
            stats = stats.split(',')
            if 'genecount' in stats:
                description.append(('genecount',"number","# genes"))
                data = self._getGeneCountHistogramData(chr)
            else:
                group = self.genomestats_file['GeneStats']
                stats_list = group['headers'][:]
                chr_numbers = ['chr1','chr2','chr3','chr4','chr5']
                chr_num = chr_numbers.index(chr.lower())
                stats_ix = [numpy.where(stats_list ==stat)[0][0] for stat in stats]
                stats_ix.sort()
                stats = stats_list[stats_ix,0]
                chr_region = group.attrs['chr_regions'][chr_num]
                stats_values = group['stats'][chr_region[0]:chr_region[1],stats_ix]
                positions = group['positions'][chr_region[0]:chr_region[1],1]

                for stat in stats:
                    stat_label = stat
                    stat_data = next((genome_stat for genome_stat in self.genome_wide_stats if genome_stat['name'] == stat), None)
                    if stat_data is not None and 'label' in stat_data:
                        stat_label = stat_data['label']
                    description.append((stat,"number",stat_label))
                stats_values_trans = stats_values.transpose().tolist()
                stats_values_trans_filtered = [map(lambda value: None if  math.isnan(value) else value,values) for values in stats_values_trans]
                data = zip(positions.tolist(),*(stats_values_trans_filtered))
            data_table = gviz_api.DataTable(description)
            data_table.LoadData(data)
            retval = {'status': 'OK','data':data_table.ToJSon()}
        except Exception,err:
            retval =  {"status":"ERROR","statustext":"%s" %str(err)}
        return retval