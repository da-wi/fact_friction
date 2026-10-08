rebuild_ctfit_from_manual <- function(ctfit, optres,
                                      savescores = TRUE,
                                      savesubjectmatrices = FALSE) {
  ctm            <- ctfit$ctstanmodel 
  #sm             <- ctfit$stanmodel
  standata       <- ctfit$standata
  stanmodeltext  <- ctfit$stanmodeltext
  data_out       <- ctfit$data
  ctdatastruct   <- ctfit$ctdatastruct
  setup <- list( 
    recompile   = ctm$recompile %||% FALSE,
    idmap       = standata$idmap,
    matsetup    = ctm$modelmats$matsetup,
    matvalues   = ctm$modelmats$matvalues,
    popsetup    = subset(ctm$modelmats$matsetup,
                         when %in% c(0L, -1L) & param > 0L),
    popvalues   = subset(ctm$modelmats$matvalues,
                         ctm$modelmats$matsetup$when %in% c(0L, -1L) &
                           ctm$modelmats$matsetup$param > 0L),
    extratforms = ctm$modelmats$extratforms
  )
  
  # 2) Make a ctStanFit$args list that looks like a normal optimize=TRUE run
  args <- ctfit$args
  args$optimize <- TRUE
  args$fit      <- TRUE

  # 3) Ensure the stanoptimis result is in the shape ctsem expects
  #    (stanoptimis already returns these fields)
  stanfit <- optres
  
  # 4) Build transformedparsfull (what summary/plots use)
  tpfull <- suppressMessages(ctsem:::stan_constrainsamples(
    sm = sm,
    standata = standata,
    savesubjectmatrices = isTRUE(savesubjectmatrices),
    samples = matrix(stanfit$rawest, 1),
    cores = 1,
    savescores = isTRUE(savescores),
    pcovn = 5000
  ))
  stanfit$transformedparsfull <- tpfull
  
  # 5) Compute Kalman outputs at the point estimate (like ctStanFit does)
  dummy_out <- list( # minimal object for ctStanKalman to read
    args = args,
    setup = setup,
    stanmodeltext = stanmodeltext,
    data = data_out,
    ctdatastruct = ctdatastruct,
    standata = standata,
    ctstanmodelbase = ctfit$ctstanmodelbase,
    ctstanmodel = ctm,
    stanmodel = sm,
    stanfit = stanfit
  )
  class(dummy_out) <- "ctStanFit"
  stanfit$kalman <- suppressMessages(ctStanKalman(dummy_out, pointest = TRUE))
  
  # 6) Final ctStanFit object
  out <- list(
    args = args,
    setup = setup,
    stanmodeltext = stanmodeltext,
    data = data_out,
    ctdatastruct = ctdatastruct,
    standata = standata,
    ctstanmodelbase = ctfit$ctstanmodelbase,
    ctstanmodel = ctm,
    stanmodel = sm,
    stanfit = stanfit
  )
  class(out) <- "ctStanFit"
  out
}