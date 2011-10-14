package com.gmi.rnaseqwebapp.client.command;

import com.gmi.rnaseqwebapp.client.dispatch.AbstractRequestBuilderCacheClientActionHandler;
import com.google.gwt.event.shared.EventBus;
import com.google.inject.Inject;
import com.gwtplatform.dispatch.client.actionhandler.caching.Cache;

public class GetPhenotypeMotionChartDataActionHandler extends AbstractRequestBuilderCacheClientActionHandler<GetPhenotypeMotionChartDataAction, GetPhenotypeMotionChartDataActionResult> {

	@Inject
	protected GetPhenotypeMotionChartDataActionHandler(Cache cache, EventBus eventBus) {
		super(GetPhenotypeMotionChartDataAction.class, cache, eventBus, true, true, false);
	}
}