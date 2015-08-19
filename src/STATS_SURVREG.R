#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "SPSS, JKP"
# version__ = "1.0.1"

# History
# 13-May-2013 Original Version


helptext="STATS SURVREG CENSORTYPE = LEFT or RIGHT* or INTERVAL or COUNTING
    EVENT=event variable name
    TIME1 = censor time variable
    TIME2 = second censor time variable for interval censoring
    INDEP = list of independent variables
    STRATA = list of stratification variables
    
/OPTIONS DISTRIBUTION = WEIBULL* or EXPONENTIAL or GAUSSIAN or 
    LOGISTIC or LOGNORMAL or LOGLOGISTIC or T
    DISTPARAM is only used for the t distribution
    SCALE = fixed scale value
    ROBUST = YES or NO*
    MISSING = LISTWISE * or INCLUDE or FAIL
    MAXITER=number PCONVERGE=relative tolerance SINGULAR = number
/OUTPUT SURVSPECS = YES or NO*
* indicates default.
CENSORTYPE, EVENT, TIME1, and INDEP are required.

Example: 
stats survreg censortype=right event=event time1=time
  indep = age quant
  /options distribution = gaussian
  /save caseresults=residuals residuals=dfbeta
  /output survspecs=yes.

CENSORTYPE specifies the type of censoring.  Interval censoring
  requires both time variables.

EVENT specifies the event variable.  The values are interpreted as
  0 = no event (alive), 1 = event occurred (dead).  
  For interval censoring, the values are 0 = right censored,
  1 = event at time, 2 = left censored, 3 = interval censored.

TIME1.  For right-censored data, the values are the follow up (endpoint) time.
  For interval censoring, the values are the interval start times.  
  
TIME2. This is only used with interval or counting process censoring.

INDEP specifies the independent variables.  Variables with a
  categorical measurement level - nominal or ordinal, will be automatically 
  converted to factors.
  
STRATA specifies stratification variables.  The same variable
    cannot appear in both the independents and strata fields.
    
DISTRIBUTION specifies the distribution of the dependent variable.
    Some distributions may be incompatible with the data, for example
    distributions with implied log transformations vs zeros

DISTPARAM  is only used for the t distribution and must be >= 3.
    
SCALE specifies a fixed scale value - any positive number or
    a zero value to indicate that it should be estimated.
    Fixed values cannot be used with multiple strata.
    
ROBUST indicates whether or not to compute robust, \"sandwich\"
    standard errors.  ROBUST is not available in some versions
    of the survival package due to an error in that package and
    will be ignored with a warning in this case.

MISSING indicates the missing value treatment.  Cases with user missing
    values can be excluded (LISTWISE), included as valid (INCLUDE), or
    the procedure can be stopped if any user missing are encountered (FAIL).
    System missing values always cause the case to be omitted.
    Missing values in the event or time variables are not allowed regardless
    of this setting.
    
MAXITER, PCONVERGE, and SINGULAR specifications control the estimation
    process.  PCONVERGE is the relative tolerance (change) below which
    convergence is determined.
    
STATS SURVREG /HELP.  prints this information and does nothing else.

This extension command requires the R programmability plug-in and 
the R eRm package.
"

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_SURVREG"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_SURVREG"))
}

    
dosurvreg = function(censortype, event, time1, time2=NULL, indep, 
    strata=NULL, distribution, 
    id=NULL, survspecs=FALSE, distparam=NULL, scale=-1, robust=FALSE,
    caseresults=NULL, resids=NULL, pred=NULL, quants=NULL,
    maxiter=30, pconverge=1e-09, singular=1e-10, missingaction="listwise"
    ) {
    #estimate survival models
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
	
	domain<-"STATS_SURVREG"
	setuplocalization(domain)
	
    procname=gtxt("Survival Regression")
    warningsprocname = gtxt("Survival Regression: Warnings")
    omsid="STATSSURVREG"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(survival), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.","survival"),dostop=TRUE)
        }
    )
    if (censortype == "interval" && is.null(time2)) {
        warns$warn(gtxt("Two time variables are required for interval censoring"), dostop=TRUE)
    }
    if (distribution == "t" && is.null(distparam)) {
        warns$warn(gtxt("The t distribution requires a degrees of freedom parameter"), dostop=TRUE)
    }
    if (distribution != "t") {
        distparam=NULL
    }
    if (!is.null(caseresults)) {
        if (is.null(c(resids, pred))) {
            warns$warn(gtxt("An output dataset was specified without any residual or prediction specifications"), 
                dostop=TRUE)
        }
        if ("quantile" %in% pred && is.null(quants)) {
            warns$warn(gtxt("Prediction quantile was specified without any quantile values"), dostop=TRUE)
        }
    } else {
        if (!is.null(c(resids, pred))) {
            warns$warn(gtxt("Residuals or predicted values were specified but no dataset name was given"),
                dostop=TRUE)
        }
    }
    if (robust && packageVersion("survival") <"2.37-5") {
        robust=FALSE
        warns$warn(gtxtf(
            "The robust option is not available in survival package version %s.  It has been set to NO.", 
            packageVersion("survival"), 80), dostop=FALSE)
    }
    # ensure that output dataset names are not in use

    alldsspecs = c(caseresults)
    if (!is.null(alldsspecs)) {
        # this api might throw an error for pending transformations
        alldatasets = spssdata.GetDataSetList()
        if (length(intersect(alldsspecs, alldatasets) > 0)) {
            warns$warn(gtxt("One or more specified output dataset names are already in use"), dostop=TRUE)
        }
    }

    if (length(union(indep, strata)) != length(indep) + length(strata)) {
        warns$warn(gtxt("The same variable cannot appear in both the independents and strata lists"),
            dostop=TRUE)
    }
    alldata = c(event, time1, time2, indep, strata)
    # if no id variable was specified, the row.label argument will be NULL and is ignored
    dta = spssdata.GetDataFromSPSS(alldata, row.label=id, missingValueToNA=TRUE, factorMode="labels")

    if (any(is.factor(dta[[time1]]), (!is.null(time2) && is.factor(dta[[time2]])), is.factor(dta[[event]]))) {
        warns$warn(gtxt("The event and time variables must have a scale measurement level"), dostop=TRUE)
    }
    for (f in strata) {
        if (!is.factor(dta[[f]])) {
            warns$warn(gtxtf("Strata variables must be categorical: %s", f), dostop=TRUE)
        }
    }
    naaction = list(listwise=na.omit, include=na.pass, fail=na.fail)[[missingaction]]
    if (is.null(time2)) {
        sur = tryCatch(Surv(dta[[time1]], dta[[event]], type=censortype),
            error=function(e) {warns$warn(e$message, dostop=TRUE)
            }
        )
    } else {
        sur = tryCatch(Surv(dta[[time1]], dta[[time2]], event=dta[[event]], type=censortype),
            error=function(e) {warns$warn(e$message, dostop=TRUE)
            }
        )
    }
    indeps = paste(indep, collapse="+")
    # build strata specs, if any, and merge with independents
    slist = list()
    for (si in length(strata)) {
        slist[si] = sprintf("strata(%s)", strata[[si]])
    }
    if (is.null(strata)) {
        rhs = indeps
    } else {
        rhs = paste(indeps, paste(slist, collapse="+"), sep="+")
    }
    fm = sprintf("sur ~ %s", rhs)

    ctrl = survreg.control(maxiter=maxiter, rel.tolerance=pconverge, toler.chol=singular)
    res = tryCatch(survreg(as.formula(fm), data=dta, na.action=naaction,
        dist=distribution, parms=distparam, scale=scale, robust=robust, control=ctrl),
        error = function(e) {
            warns$warn(e$message, dostop=FALSE)
            return(NULL)
        }
    )
    
    # if failed and survival interpretation requested, display that much and stop
    if (is.null(res)) {
        if (survspecs) {
            StartProcedure(gtxt(procname), omsid)
            survtable(sur, dta, event, time1, time2)
        }
        warns$warn(dostop=TRUE)
    }
            
    ressum = summary(res)
    StartProcedure(gtxt("Survival Regression"), "STATSSURVREG")
    summarylabels=list(
        gtxt("Event Variable"),
        gtxt("Censoring Type"),
        gtxt("Time"),
        gtxt("Time2"),
        gtxt("Distribution"),
        gtxt("Robust Estimation"),
        gtxt("Output Dataset"),
        gtxt("Missing Values Treatment"),
        gtxt("N"),
        gtxt("Number of Iterations"),
        gtxt("Scale"),
        gtxt("Log Likelihood (Model)"),
        gtxt("Log Likelihood (Intercept Only)"),
        gtxt("Chi-Squared"),
        gtxt("D.F."),
        gtxt("P Value")
    )
    summaryvalues = list(
        event,
        censortype,
        time1,
        ifelse(is.null(time2), gtxt("--None--"), time2),
        ifelse(is.null(distparam), distribution, gtxtf("%s(%s)", distribution, distparam)),
        ifelse(ressum$robust, gtxt("Yes"), gtxt("No")),
        ifelse(is.null(caseresults), gtxt("--None--"), caseresults),
        missingaction,
        ressum$n,
        ressum$iter,
        ifelse(scale <= 0, res$scale, gtxtf("Fixed at %s", scale)),
        res$loglik[[2]],
        res$loglik[[1]],
        ressum$chi,
        ressum$df - ressum$idf,
        1 - pchisq(ressum$chi, ressum$df - ressum$idf)
    )
    names(summaryvalues) = summarylabels
    summarydf = data.frame(cbind(summaryvalues))
    colnames(summarydf) = gtxt("Values")
    spsspivottable.Display(summarydf, title=gtxt("Survival Regression Summary"), 
        templateName="SURVREGSUMMARY",
        caption=gtxt("Results computed by R survival package function survreg"),
        isSplit=FALSE,
        format=formatSpec.Coefficient
    )
    
    resdf = data.frame(ressum$table)
    names(resdf) = c(gtxt("Value"), gtxt("Std. Error"), gtxt("Z"), gtxt("P"))
    caption = gtxtf("Event variable: %s", event)
    if (!is.null(strata)) {
        caption = sprintf(gtxt("%s.  Strata variables: %s"), caption, paste(strata, sep=", "))
    }
    spsspivottable.Display(resdf, title=gtxt("Survival Regression"), templateName="SURVREGRES",
        caption=caption, isSplit=FALSE
    )
    
    if (survspecs) {
        survtable(sur, dta, event, time1, time2)
    }
    spsspkg.EndProcedure()
    
    if (!is.null(caseresults)) {
        docaseresults(caseresults, res, resids, pred, quants, row.names(dta), id, warns)
    }
    warns$warn()  # ensure warnings are displayed
}
survtable = function(sur, dta, event, time1, time2) {
    # Display partial survival interpretation
      # The Surv object actually has two columns and gets special formatting
      # when displayed, but trickery is required for the pivot table formatting
      limit = min(nrow(sur), 50)
      surdf = data.frame(sur[1:limit,])
      attr(surdf, "names") = gtxt("Survival")
      surdf = sapply(surdf, as.character)
      for (v in c(event, time1, time2)) {
          surdf = data.frame(surdf, dta[1:limit, v])
          names(surdf)[length(surdf)] = v
      }
      spsspivottable.Display(surdf, title=gtxtf("Survival Values (First %s Cases)", limit),
          templateName="SURVVALUES",
          isSplit=FALSE
      )
}        
        
labels = list(res_response=gtxt("Residuals - Response"),
    res_deviance = gtxt("Residuals - Deviance"),
    res_working = gtxt("Residuals - Working"),
    res_ldcase = gtxt("Residuals - Ldcase"),
    res_ldresp = gtxt("Residuals - Ldresponse"),
    res_ldshape = gtxt("Residuals - Ldshape"),
    res_dfbeta = gtxt("Dfbeta"),
    res_dfbetas = gtxt("Stdized Dfbeta"),
    pred_response = gtxt("Predicted - Response"),
    pred_linear = gtxt("Predicted - Linear"),
    pred_quantile = gtxt("Predicted - Quantile"),
    pred_uquantile = gtxt("Predicted - Uquantile")
)

docaseresults = function(caseresults, res, resids, pred, quants, id, idname, warns) {

    dict = list()
    if (is.null(idname) || length(id) != nrow(res$y)) {
        df = data.frame(row.names=1:length(res$y))  # empty data frame with the right number of rows
        dict["ID"] = list(c("ID", "", 0, "F8.0", "nominal"))
    } else {
        df = data.frame(row.names=id)  # empty data frame with the right number of rows
        dict[["ID"]] = spssdictionary.GetDictionaryFromSPSS(idname)
        dict[["ID"]][1,] = "ID"     # ensure that we can't have a duplicate name
    }
    resultcount = 0

    # first build residual results; then prediction results
    for (rtype in resids) {
        e = tryCatch(residuals(res, type=rtype), # some types just fail with this function for some models :-)
                error = function(ea) {
                    warns$warn(gtxtf("Residual statistic %s not available for this model", rtype), dostop=FALSE)
                    return(NULL)
                }
            )
        # add successfull results to the data frame

        if (!is.null(e)) {
            if (!(rtype %in% list("dfbeta", "dfbetas"))) {  # dfbeta/s are multicolumn
                resultcount = resultcount + 1
                df = data.frame(df, e)
                itemname = paste("res", rtype, sep ="_")
                dict[itemname] = list(c(itemname, labels[[itemname]], 0, "F8.2", "scale"))
            } else {
                for (v in 1:ncol(e)) {
                    resultcount = resultcount + 1
                    df = data.frame(df, e[,v])
                    itemname = paste(rtype, v, sep="_")  # improve this
                    if (v < ncol(e)) {
                        itemlabel = names(res$coefficients)[[v]]
                    } else {
                        itemlabel = gtxt("Log(scale)")
                    } 
                    dict[itemname] = list(c(itemname, itemlabel, 0, "F8.2", "scale"))
                }
            }
        }
    }
    
    # predictions
    for (rtype in pred) {
        if (!(rtype %in% list("quantile", "uquantile"))) {   # quantile output can be multicolumn
            e = tryCatch(predict(res, type=rtype, se.fit=TRUE), 
                    error = function(ea) {
                        warns$warn(gtxtf("Prediction statistic %s not available for this model", rtype), dostop=FALSE)
                        return(NULL)
                    }
                )
            if (!is.null(e)){
                resultcount = resultcount + 1
                df = data.frame(df, e)
                itemname = paste("pred", rtype, sep="_")
                dict[itemname] = list(c(itemname, labels[[itemname]], 0, "F8.2", "scale"))
                resultcount = resultcount + 1   # std error
                itemname = paste(itemname, "SE", sep="_")
                dict[itemname] = list(c(itemname, itemname, 0, "F8.2", "scale"))
            }
        } else {  # quantile and uquantile.  se.fit does not work for these (undocumented)
            e = tryCatch(predict(res, type = rtype, se.fit=FALSE, p=quants),
                    error = function(ea) {
                        warns$warn(gtxtf("Prediction statistic %s not available for this model", rtype), 
                            dostop=FALSE)
                        return(NULL)
                    }
                )
            if (!is.null(e)) {
                df = data.frame(df, e)   # as many columns as quantile values
                for (q in quants) {
                    resultcount = resultcount + 1
                    itemname = paste("pred", 
                        ifelse(rtype == "quantile", gtxt("quantile"), gtxt("uquantile")),
                        q,
                        sep = "_"
                        )
                    dict[itemname] = list(c(itemname, itemname, 0, "F8.2", "scale"))
                }
            }
        }
    }
    
    # finally! build the dataset
    if (resultcount == 0) {
        warns$warn(gtxt("No casewise results could be computed"), dostop=TRUE)
    }

    tryCatch({
        spssdictionary.SetDictionaryToSPSS(caseresults, 
            spssdictionary.CreateSPSSDictionary(dict))
        spssdata.SetDataToSPSS(caseresults, data.frame(row.names(df), df))
        },
        error=function(e) {warn.warn(e)
          warns$warn(gtxt("Failed to create caseresults dataset."), dostop=FALSE)}
    )
    spssdictionary.EndDataStep()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spss.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 


Run = function(args) {
    #Execute the STATS SURVREG command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("CENSORTYPE", subc="", ktype="str", var="censortype",
            vallist=list("right", "left", "interval", "counting")),
        spsspkg.Template("EVENT", subc="",  ktype="existingvarlist", var="event",
            islist=FALSE),
        spsspkg.Template("TIME1", subc="", ktype="existingvarlist", var="time1",
            islist=FALSE),
        spsspkg.Template("TIME2", subc="",  ktype="existingvarlist", var="time2"),
        spsspkg.Template("INDEP", subc="", ktype="existingvarlist", var="indep", islist=TRUE),
        spsspkg.Template("STRATA", subc="", ktype="existingvarlist", var="strata", islist=TRUE),
        spsspkg.Template("ID", subc="", ktype="existingvarlist", var="id"),
        
        spsspkg.Template("DISTRIBUTION", subc="OPTIONS", ktype="str", var="distribution",
            vallist=list("weibull", "exponential", "gaussian", "logistic",
            "lognormal", "loglogistic", "extreme", "loggaussian", "rayleigh", "t")),
        spsspkg.Template("DISTPARAM", subc="OPTIONS", ktype="float", var="distparam"),
        spsspkg.Template("SCALE", subc="OPTIONS", ktype="float", var="scale"),
        spsspkg.Template("ROBUST", subc="OPTIONS", ktype="bool", var="robust"),
        spsspkg.Template("MAXITER", subc="OPTIONS", ktype="int", var="maxiter",
            vallist=list(1)),
        spsspkg.Template("PCONVERGE", subc="OPTIONS", ktype="float", var="pconverge",
            vallist=list(0)),
        spsspkg.Template("SINGULAR", subc="OPTIONS", ktype="float", var="singular",
            vallist=list(0)),
        spsspkg.Template("MISSING", subc="OPTIONS", ktype="str", var="missingaction",
            vallist=list("listwise", "fail")),
        
        spsspkg.Template("SURVSPECS", subc="OUTPUT", ktype="bool", var="survspecs"),
        
        spsspkg.Template("CASERESULTS", subc="SAVE", ktype="varname", var="caseresults"),
        #,"ldcase","ldresp","shape" sometimes do not work in survreg
        spsspkg.Template("RESIDUALS", subc="SAVE", ktype="str", var="resids", islist=TRUE,
            vallist = list("response","deviance","working", "dfbeta","dfbetas","ldcase","ldresp","ldshape")),  
        spsspkg.Template("PREDICT", subc="SAVE", ktype="str", var="pred", islist=TRUE,
            vallist=list("response","linear","quantile","uquantile")),
        spsspkg.Template("QUANTILES", subc="SAVE",ktype="float", var="quants", islist=TRUE,
            vallist=list(0.0001, .9999))
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "dosurvreg")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}