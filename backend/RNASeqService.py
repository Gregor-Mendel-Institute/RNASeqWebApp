'''
Created on Sep 12, 2011

@author: uemit.seren
'''
import gviz_api, os
import tables
from datetime import datetime
from JBrowseDataSource import DataSource as JBrowseDataSource 
import cherrypy
import cPickle
from variation.src.rnaseq_records import RNASeqRecords 

class RNASeqService:
    base_path = "/net/gmi.oeaw.ac.at/gwasapp/rnaseq-web/"
    base_path_jbrowse = base_path
    base_path_datasets = base_path + "datasets/"
    base_jbrowse_path = base_path_jbrowse + "jbrowse_1.2.1/"
    track_folder = "TAIR10"
    
    #tracks/Chr%s/TAIR10/"
    __datasource = None
    _lazyArrayChunks = [{}, {}, {}, {}, {}]
    hdf5_filename = base_path + "rnaseq_data.hdf5"
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
        self.phenotypes = cPickle.load(phenotype_file)
        phenotype_file.close()
            
    
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
    
        
    
    def _getPhenotypes(self, range_start,range_length,name='', chr='', start='', end=''):
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
        retval['phenotype'] = phenotype;
        retval['histogramdataTable'] = self._getCombinedHistogramDataTable(id)
        return retval
    
    @cherrypy.expose
    @cherrypy.tools.json_out()
    def getPhenotypeMotionChartData(self,phenotype,environment,dataset,transformation):
        return {"data":self._getPhenotypeExplorerData(phenotype,environment,dataset,transformation)}

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
                data = zip(association_result[i]['position'],association_result[i]['score'])
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