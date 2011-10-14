package com.gmi.rnaseqwebapp.client.command;

import java.util.List;

import com.gmi.rnaseqwebapp.client.dto.Cofactor;
import com.gmi.rnaseqwebapp.client.dto.ResultData;
import com.gwtplatform.dispatch.shared.Result;

public class GetGWASDataActionResult implements Result {
	private final ResultData resultData;
	
	public GetGWASDataActionResult(final ResultData resultData) {
		this.resultData = resultData;
	}
	
	public ResultData getResultData() {
		return this.resultData;
	}


}
