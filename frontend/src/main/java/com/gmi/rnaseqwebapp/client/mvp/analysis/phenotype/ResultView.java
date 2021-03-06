package com.gmi.rnaseqwebapp.client.mvp.analysis.phenotype;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import at.gmi.nordborglab.widgets.geneviewer.client.datasource.DataSource;
import at.gmi.nordborglab.widgets.geneviewer.client.datasource.Gene;
import at.gmi.nordborglab.widgets.geneviewer.client.datasource.impl.GeneSuggestion;
import at.gmi.nordborglab.widgets.geneviewer.client.datasource.impl.ServerSuggestOracle;
import at.gmi.nordborglab.widgets.gwasgeneviewer.client.GWASGeneViewer;

import com.gmi.rnaseqwebapp.client.dto.Cofactor;
import com.gmi.rnaseqwebapp.client.dto.Phenotype;
import com.gmi.rnaseqwebapp.client.resources.CellTableResources;
import com.gmi.rnaseqwebapp.client.ui.ResizeableFlowPanel;
import com.google.gwt.dom.client.AnchorElement;
import com.google.gwt.event.dom.client.ChangeEvent;
import com.google.gwt.event.logical.shared.SelectionEvent;
import com.google.gwt.event.logical.shared.SelectionHandler;
import com.google.gwt.event.logical.shared.ValueChangeEvent;
import com.google.gwt.uibinder.client.UiBinder;
import com.google.gwt.uibinder.client.UiField;
import com.google.gwt.uibinder.client.UiHandler;
import com.google.gwt.user.cellview.client.CellTable;
import com.google.gwt.user.cellview.client.TextColumn;
import com.google.gwt.user.client.ui.CheckBox;
import com.google.gwt.user.client.ui.HTMLPanel;
import com.google.gwt.user.client.ui.IsWidget;
import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.SuggestBox;
import com.google.gwt.user.client.ui.SuggestBox.DefaultSuggestionDisplay;
import com.google.gwt.user.client.ui.SuggestOracle;
import com.google.gwt.user.client.ui.SuggestOracle.Suggestion;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.view.client.HasData;
import com.google.gwt.visualization.client.DataTable;
import com.google.gwt.visualization.client.DataView;
import com.google.gwt.visualization.client.visualizations.corechart.LineChart;
import com.google.gwt.visualization.client.visualizations.corechart.Options;
import com.google.inject.Inject;
import com.gwtplatform.mvp.client.ViewWithUiHandlers;

public class ResultView extends ViewWithUiHandlers<ResultUiHandlers> implements ResultPresenter.MyView {

	private final Widget widget;

	public interface Binder extends UiBinder<Widget, ResultView> {
	}
	
	@UiField AnchorElement download_link; 
	@UiField ResizeableFlowPanel gwas_container;
	@UiField(provided = true) CellTable<Cofactor> cofactorTable;
	@UiField(provided=true) LineChart statistic_chart;
	@UiField ListBox statistic_type;
	@UiField HTMLPanel statistic_container;
	@UiField HTMLPanel model_container;
	@UiField CheckBox full_check;
	@UiField CheckBox genetic_check;
	@UiField CheckBox environ_check;
	private final DataSource geneDataSource;
	private final CellTableResources cellTableResources;
	private String[] colors = {"blue", "green", "red", "cyan", "purple"};
	private String[] stacked_colors = {"rgba(0,0,255,0.4)","rgba(255,0,0,0.4)","rgba(0,255,0,0.4)"};
	private String[] gene_mark_colors = {"red", "red", "blue", "red", "green"};
	protected List<GWASGeneViewer> gwasGeneViewers = new ArrayList<GWASGeneViewer>();
	protected Phenotype phenotype;
	@UiField(provided=true)	final SuggestBox searchGene;

	@Inject
	public ResultView(final Binder binder, final DataSource geneDataSource, final CellTableResources cellTableResources) {
		this.geneDataSource = geneDataSource;
		this.cellTableResources = cellTableResources;
		searchGene = new SuggestBox(new ServerSuggestOracle(geneDataSource,5));
		searchGene.getElement().setAttribute("placeHolder", "Search gene");
		((DefaultSuggestionDisplay)searchGene.getSuggestionDisplay()).setAnimationEnabled(true);
		searchGene.addSelectionHandler(new SelectionHandler<SuggestOracle.Suggestion>() {
			
			@Override
			public void onSelection(SelectionEvent<Suggestion> event) {
				GeneSuggestion suggestion =  (GeneSuggestion)event.getSelectedItem();
				Gene gene = suggestion.getGene();
				Gene phenotypeGene = getGeneFromPhenotype();
				GWASGeneViewer viewer = getGWASGeneViewer(gene.getChromosome());
				if (viewer != null)
				{
					viewer.clearDisplayGenes();
					if (viewer.getChromosome().equals(phenotypeGene.getChromosome())) {
						viewer.addDisplayGene(getGeneFromPhenotype());
					}
					viewer.addDisplayGene(gene);
					viewer.refresh();
				}
			}
		});
		
		cofactorTable = new CellTable<Cofactor>(15,cellTableResources);
		cofactorTable.addColumn(new TextColumn<Cofactor>() {

			@Override
			public String getValue(Cofactor object) {
				return object.getStep().toString();
			}
		},"Step");
		cofactorTable.addColumn(new TextColumn<Cofactor>() {

			@Override
			public String getValue(Cofactor object) {
				if (object.getChr() != null)
					return object.getChr().toString();
				else return "";
			}
		},"Chr.");
		cofactorTable.addColumn(new TextColumn<Cofactor>(){

			@Override
			public String getValue(Cofactor object) {
				if (object.getPos() != null)
					return object.getPos().toString();
				else
					return "";
			}
			
		}, "Position");
		cofactorTable.addColumn(new TextColumn<Cofactor>(){

			@Override
			public String getValue(Cofactor object) {
				return Double.toString(object.getBic());
			}
			
		}, "BIC");
		cofactorTable.addColumn(new TextColumn<Cofactor>(){

			@Override
			public String getValue(Cofactor object) {
				return Double.toString(object.getEbic());
			}
			
		}, "e-BIC");
		cofactorTable.addColumn(new TextColumn<Cofactor>(){

			@Override
			public String getValue(Cofactor object) {
				return Double.toString(object.getMbic());
			}
			
		}, "m-BIC");
		cofactorTable.addColumn(new TextColumn<Cofactor>(){

			@Override
			public String getValue(Cofactor object) {
				return Double.toString(object.getMbonf());
			}
			
		}, "mbonf");
		cofactorTable.addColumn(new TextColumn<Cofactor>(){

			@Override
			public String getValue(Cofactor object) {
				return Double.toString(object.getPseudoHeritability());
			}
			
		}, "pseudo-heritability");
		
		statistic_chart = new LineChart(DataTable.create(), createStatsticChartOptions());
		widget = binder.createAndBindUi(this);
		
		if (statistic_type.getItemCount() == 0)
		{
			statistic_type.addItem("", "");
			statistic_type.addItem("BIC", "1");
			statistic_type.addItem("mBIC", "2");
			statistic_type.addItem("eBIC", "3");
			statistic_type.addItem("mBonf", "4");
			statistic_type.addItem("pseudo-heritability", "5");
		}
	}

	@Override
	public Widget asWidget() {
		return widget;
	}

	@Override
	public void setDownloadURL(String url) {
		download_link.setHref(url);
	}

	@Override
	public void drawAssociationCharts(List<DataTable> dataTables,
			List<Cofactor> cofactors, List<Integer> chrLengths,
			double maxScore, double bonferroniThreshold,boolean isStacked,
			Phenotype phenotype) {
		clearAssociationCharts();
		this.phenotype = phenotype;
		Gene gene = getGeneFromPhenotype();
		Integer i = 1;
		java.util.Iterator<DataTable> iterator = dataTables.iterator();
		List<Cofactor> cofactors_to_add = new ArrayList<Cofactor>();
		if (cofactors.size() > 0)
			cofactors_to_add.addAll(cofactors.subList(1, cofactors.size()));
		String[] color = null;
		if (isStacked) {
			color = stacked_colors;
			color = new String[] {colors[0],colors[1],colors[2]};
		}
		while(iterator.hasNext())
		{
			GWASGeneViewer chart =null;
			DataTable dataTable = iterator.next();
			if (!isStacked)
				color = new String[] {colors[i%colors.length]};
			String gene_marker_color = gene_mark_colors[i%gene_mark_colors.length];
			if (gwasGeneViewers.size() >= i)
				chart = gwasGeneViewers.get((i-1));
			if (chart == null)
			{
				chart = new GWASGeneViewer("Chr"+i.toString(), color, gene_marker_color,geneDataSource);
				chart.setGeneInfoUrl("http://arabidopsis.org/servlets/TairObject?name={0}&type=gene");
				gwasGeneViewers.add(chart);
				gwas_container.add((IsWidget)chart);
			}
			else {
				chart.clearDisplayGenes();
				chart.setColor(color);
			}
			chart.clearSelection();
			Iterator<Cofactor> itr = cofactors_to_add.iterator();
			while(itr.hasNext()){
				Cofactor cofactor = itr.next();
				if (cofactor.getChr() == i)
				{
					chart.addSelection(GWASGeneViewer.getSelectionFromPos(dataTable, cofactor.getPos()));
					itr.remove();
				}
			}
			if (chart.getChromosome().equals(gene.getChromosome()))
				chart.addDisplayGene(gene);
			chart.draw(dataTable,maxScore,0,chrLengths.get(i-1),bonferroniThreshold);
			i++;
		}
		
	}
	
	protected GWASGeneViewer getGWASGeneViewer(String chromosome) {
		for (GWASGeneViewer gwasGeneViewer : gwasGeneViewers) {
			
			if (gwasGeneViewer.getChromosome().equals(chromosome)) 
				return gwasGeneViewer;
		}
		return null;
	}


	public void clearAssociationCharts() {
		/*for (GWASGeneViewer gwasgeneViewer: gwasGeneViewers) {
			gwasgeneViewer.destroy();
		}
		gwas_container.clear();
		gwasGeneViewers.clear();*/
	}
	
	private Options createStatsticChartOptions() {
		Options options = Options.create();
	    options.setWidth(600);
	    options.setHeight(400);
	    return options;
	}
	
	@Override
	public void drawStatisticPlots(DataView view) {
		if (statistic_chart == null) {
			statistic_chart = new LineChart(view, createStatsticChartOptions());
			statistic_container.add(statistic_chart);
		}
		else
			statistic_chart.draw(view, createStatsticChartOptions());
	}

	@Override
	public HasData<Cofactor> getCofactorDisplay() {
		return cofactorTable;
	}
	
	@UiHandler("statistic_type")
	public void onChangeStatisticType(ChangeEvent event) {
		if (statistic_type.getSelectedIndex() <=0 )
			statistic_chart.setVisible(false);
		else
			getUiHandlers().loadStatisticChart(statistic_type.getValue(statistic_type.getSelectedIndex()));
	}

	@Override
	public void reset() {
		detachCharts();
		clearAssociationCharts();
	}

	@Override
	public void detachCharts() {
		statistic_type.setSelectedIndex(0);
		statistic_container.clear();
		statistic_chart = null;
	}
	
	@UiHandler({"full_check","genetic_check","environ_check"})
	public void onClickCheckBox(ValueChangeEvent<Boolean> event) 
	{
		CheckBox checkbox = (CheckBox)event.getSource();
		int id = 0;
		if (checkbox == full_check) {
			id = 0;
		}
		else if (checkbox == genetic_check) {
			id = 1;
		}
		else if (checkbox == environ_check) {
			id = 2;
		}
		for (GWASGeneViewer gwasgeneViewer: gwasGeneViewers) {
			gwasgeneViewer.setScatterChartVisibilityForSeries(id, event.getValue());
		}
	}

	@Override
	public void setIsStacked(boolean isStacked) {
		model_container.setVisible(isStacked);
		full_check.setValue(true);
		genetic_check.setValue(true);
		environ_check.setValue(true);
	}
	
	protected Gene getGeneFromPhenotype() {
		if (phenotype != null) {
			return new Gene(phenotype.getName(),"Chr"+phenotype.getChr().toString(),phenotype.getStart(),phenotype.getEnd(),"");
		}
		return null;
	}
}