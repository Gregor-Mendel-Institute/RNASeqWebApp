package com.gmi.rnaseqwebapp.client.command;


import com.gmi.rnaseqwebapp.client.dispatch.AbstractRequestBuilderCacheClientActionHandler;
import com.google.web.bindery.event.shared.EventBus;
import com.google.inject.Inject;
import com.gwtplatform.dispatch.client.actionhandler.caching.Cache;

public class GetGWASDataActionHandler extends AbstractRequestBuilderCacheClientActionHandler<GetGWASDataAction, GetGWASDataActionResult> {

	
	@Inject
	protected GetGWASDataActionHandler(Cache cache,EventBus eventBus) {
		super(GetGWASDataAction.class, cache,eventBus,true,true,false);
	}

}
