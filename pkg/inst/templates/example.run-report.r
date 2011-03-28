\documentclass{article}

\usepackage{geometry} 
\usepackage{graphicx}
\usepackage{longtable}
\usepackage{pdflscape}
\usepackage{fancyhdr}
\usepackage[iso]{isodateo}

\geometry{tmargin=2cm,bmargin=2.5cm,lmargin=2cm,rmargin=2cm}
\pagestyle{fancy}
\headheight 25pt

\begin{document}

\lhead{Data interpretation report}
\rhead{Acquisition Time: \Sexpr{attr(immunoassay.environment$l.run,"Date")} }
\rfoot{Page: \thepage}
\cfoot{Report Date: \today}

\setcounter{secnumdepth}{-1 }

<<results=hide, echo=FALSE>>=
	# Load libraries
	library(xtable)
	library(Cairo)

	# Prepare data for table
	tab1 = tab2 = pfile = vector("list", immunoassay.environment$n)
	for (i2 in 1:immunoassay.environment$n) {
		tab        = print(immunoassay.environment$l.run)
		tab1[[i2]] = tab[,c("Type", paste(c("MFI","CTS","conc","pred"), immunoassay.environment$l.analytes[i2], 
			sep="."))]
        tab1[[i2]]$fitted = NA
        tab1[[i2]]$fitted[tab1[[i2]]$Type=="Standard"] = immunoassay.environment$fits[[i2]]$data$use
		tab        = project.means(immunoassay.environment$l.run)
		tab2[[i2]] = tab[,c("SPL", "Type", paste(c("MFI","CTS","conc","pred"),
			immunoassay.environment$l.analytes[i2], sep="."), paste("pred", immunoassay.environment$l.analytes[i2], "cv", sep="."))]
		tab2[[i2]]$prec = (1- (tab[,paste("conc", immunoassay.environment$l.analytes[i2], sep=".")] / 
            tab[,paste("pred", immunoassay.environment$l.analytes[i2], sep=".")])) *100
		
		# Set file name
		i3   	    = gsub("_", "-", attr(immunoassay.environment$l.run, "file"), fixed=T)
		pfile[[i2]] = paste("pdf/", i3, "-", i2,".pdf", sep="")
		
		# Make the plot
		CairoPDF(file=paste(immunoassay.environment$path, "/Tex/", pfile[[i2]], sep=""), width=14, height=11)
			par(mfrow=c(2,2))
            plot(immunoassay.environment$fits[[i2]])
			plot(immunoassay.environment$fits[[i2]], type="resid", bg="#bbccaa99")
            plot(immunoassay.environment$l.run, analyte=i2, type="cv")
            plot(immunoassay.environment$l.run, analyte=i2, type="accuracy")
		dev.off()
		}
@

<<results=tex, echo=FALSE>>= 
    for (i2 in 1:immunoassay.environment$n) {
    
    # Header & footer
    cat("\\begin{center}\n", sep="")
    cat("   \\section{\\LARGE{Plate: ", i3, "}}\n", sep="")
	cat("   \\section{\\large{Analyte: ", attr(immunoassay.environment$l.run, "Analytes")[i2], "}}\n", sep="")
    cat("\\end{center}\n", sep="")
    cat("\n\\lfoot{Analyte: ", attr(immunoassay.environment$l.run, "Analytes")[i2], " }\n\n", sep="")

    # Run information
    cat("\\begin{minipage}{0.55\\textwidth}\\footnotesize\n", sep="")
	cat("  \\begin{flushleft}\n", sep="")
	cat("   \\emph{Platform: }\n    ", attr(immunoassay.environment$l.run,"Assay")[1], "\\\\ \n", sep="")
    cat("   \\emph{Instrument S/N: }\n    ", attr(immunoassay.environment$l.run,"Assay")[2], "\\\\ \n", sep="")
    cat("   \\emph{Operator: }\n    ", attr(immunoassay.environment$l.run,"Operator"), "\\\\ \n", sep="")
	cat("  \\end{flushleft}\n", sep="")
    cat("\\end{minipage}\n", sep="")
    cat("\\begin{minipage}{0.37\\textwidth}\\footnotesize\n", sep="")
	cat("  \\begin{flushright}\n", sep="")
    cat("   \\emph{Standards \\& Controls Lot \\#: }\n    ", gsub("_", "-", attr(immunoassay.environment$l.run,"Kit")[1], fixed=T), "\\\\ \n", sep="")
	cat("   \\emph{Background MFI: }\n    ", attr(immunoassay.environment$l.run,"Background")[i2], "\\\\ \n", sep="")
	cat("  \\end{flushright}\n", sep="")
    cat("\\end{minipage}\n", sep="")

    cat("\n\\vspace{0.2cm}\n\n", sep="")

    # Fit information
    cat("\\begin{minipage}{0.55\\textwidth}\\footnotesize\n", sep="")
    cat("  \\begin{flushleft}\n", sep="")
    cat("   \\emph{Fit $R^2 = $}\n    ", round(check(immunoassay.environment$fits[[i2]])$r.squared,4), "\\\\ \n", sep="")
    cat("   \\emph{$S_{y|x} = $}\n    ", round(check(immunoassay.environment$fits[[i2]])$Syx,4), "\\\\ \n", sep="")
    cat("   \\emph{Residual deviance = }\n    ", round(check(immunoassay.environment$fits[[i2]])$sigma,2), "\\\\ \n", sep="")
    cat("   \\end{flushleft}\n", sep="")
    cat("  \\end{minipage}\n", sep="")
    cat("\\begin{minipage}{0.37\\textwidth}\\footnotesize\n", sep="")
    cat("  \\begin{flushright}\n", sep="")
    cat("   \\emph{Fit type: }\n    ", paste(immunoassay.environment$fits[[i2]]$model[1], ".", 
        immunoassay.environment$fits[[i2]]$model[2], " with ", immunoassay.environment$fits[[i2]]$model[3], 
        " weights.", sep=""), "\\\\ \n", sep="")
    cat("   \\emph{Mean Calibrator Inaccuracy [\\%]: }\n    ", round(check(immunoassay.environment$fits[[i2]])$St.error[2],1), "\\\\ \n", sep="")
    cat("   \\emph{Mean QC Inaccuracy [\\%]: }\n    ", round(check(immunoassay.environment$fits[[i2]])$QC.error[2],1), "\\\\ \n", sep="")
    cat("  \\end{flushright}\n", sep="")
    cat("\\end{minipage}\n\n", sep="")

    # Add plot
    cat("\\begin{center}\n", sep="")
	cat("   \\includegraphics[width=1.01\\textwidth]{\"", substr(pfile[[i2]], 0, nchar(pfile[[i2]])-4), "\"}\n", sep="")
    cat("\\end{center}\n\n", sep="")
 
	# Print the table of means  - standards
	cat("\\begin{center}\\footnotesize \r\n Table of means: Standards \r\n \\begin{longtable}{p{1.5cm}r *{6}{c}} \r\n", sep="")
	cat("\\hline \r\n  SPL & Type & MFI & Counts & Expected & Measured & Precision & Inaccuracy \\\\ \r\n", sep="") 
	cat("\\hline \r\n   \\endhead \r\n    \\hline \r\n    \\endfoot \r\n\n 	   \\hline \\hline \r\n	   \\endlastfoot \r\n", sep="")
	print(xtable(tab2[[i2]][tab2[[i2]]$Type=="Standard",]), include.rownames=F, include.colnames=F, hline.after=NULL, only.contents=T)
	cat("\r\n \\hline \r\n", "Mean & & & & & & ", round(mean(abs(tab2[[i2]][tab2[[i2]]$Type=="Standard", 7]), na.rm=T),2), " & ",
		round(mean(abs(tab2[[i2]][tab2[[i2]]$Type=="Standard", 8]), na.rm=T),2), " \\\\"," \r\n", sep="")
	cat("\r\n	\\end{longtable} \r\n  \\end{center}", "\r\n", "\\newpage", "\r\n", sep="")
	cat("\r\n", "\r\n", sep="")
	
	# Print the table of means  - QCs
    if (nrow(tab2[[i2]][tab2[[i2]]$Type=="QC",])>0) {
        cat("\\begin{center}\\footnotesize \r\n Table of means: QCs \r\n \\begin{longtable}{p{1.5cm}r *{6}{c}} \r\n", sep="")
        cat("\\hline \r\n  SPL & Type & MFI & Counts & Expected & Measured & Precision & Inaccuracy \\\\ \r\n", sep="") 
        cat("\\hline \r\n   \\endhead \r\n    \\hline \r\n    \\endfoot \r\n\n 	   \\hline \\hline \r\n	   \\endlastfoot \r\n", sep="")
        print(xtable(tab2[[i2]][tab2[[i2]]$Type=="QC",]), include.rownames=F, include.colnames=F, hline.after=NULL, only.contents=T)
        cat("\r\n \\hline \r\n", "Mean & & & & & & ", round(mean(abs(tab2[[i2]][tab2[[i2]]$Type=="QC", 7]), na.rm=T),2), " & ",
            round(mean(abs(tab2[[i2]][tab2[[i2]]$Type=="QC", 8]), na.rm=T),2), " \\\\"," \r\n", sep="")
        cat("\r\n	\\end{longtable} \r\n  \\end{center}", "\r\n", sep="")
        cat("\r\n", "\r\n", sep="")
        }

	# Print the table of means  - Pools
	cat("\\begin{center}\\footnotesize \r\n Table of means: Pools \r\n \\begin{longtable}{p{1.5cm}r *{6}{c}} \r\n", sep="")
	cat("\\hline \r\n  SPL & Type & MFI & Counts & Expected & Measured & Precision \\\\ \r\n", sep="") 
	cat("\\hline \r\n   \\endhead \r\n    \\hline \r\n    \\endfoot \r\n\n 	   \\hline \\hline \r\n	   \\endlastfoot \r\n", sep="")
	print(xtable(tab2[[i2]][tab2[[i2]]$Type=="Pool", 1:7]), include.rownames=F, include.colnames=F, hline.after=NULL, only.contents=T)
	cat("\r\n \\hline \r\n", "Mean & & & & & & ", round(mean(abs(tab2[[i2]][tab2[[i2]]$Type=="Pool", 7]), na.rm=T),2), " \\\\"," \r\n", sep="")
	cat("\r\n	\\end{longtable} \r\n  \\end{center}", "\r\n", sep="")
	cat("\r\n", "\r\n", sep="")

	# Print the table of means  - Samples
	cat("\\begin{center}\\footnotesize \r\n Table of means: Samples \r\n \\begin{longtable}{p{2.5cm}r *{6}{c}} \r\n", sep="")
	cat("\\hline \r\n  SPL & Type & MFI & Counts & Expected & Measured & Precision \\\\ \r\n", sep="") 
	cat("\\hline \r\n   \\endhead \r\n    \\hline \r\n    \\endfoot \r\n\n 	   \\hline \\hline \r\n	   \\endlastfoot \r\n", sep="")
	print(xtable(tab2[[i2]][tab2[[i2]]$Type %in% c("Sample","Retest"), 1:7]), include.rownames=F, 
        include.colnames=F, hline.after=NULL, only.contents=T)
	cat("\r\n \\hline \r\n", "Mean & & & & & & ", round(mean(abs(tab2[[i2]][tab2[[i2]]$Type %in% c("Sample",
        "Retest"), 7]), na.rm=T),2), " \\\\"," \r\n", sep="")
	cat("\r\n	\\end{longtable} \r\n  \\end{center}", "\r\n", sep="")
	cat("\r\n", "\r\n", sep="")
	
	# Print the table of duplicates - all samples
	cat("\\begin{center}\\footnotesize \r\n Table of duplicate results \r\n \\begin{longtable}{p{4.0cm}r *{6}{c}} \r\n", sep="")
	cat("\\hline \r\n  SPL & Type & MFI & Counts & Expected & Measured & Used \\\\ \r\n", sep="") 
	cat("\\hline \r\n   \\endhead \r\n    \\hline \r\n    \\endfoot \r\n\n 	   \\hline \\hline \r\n	   \\endlastfoot \r\n", sep="")
	print(xtable(tab1[[i2]]), include.colnames=F, hline.after=NULL, only.contents=T)
	cat("\r\n	\\end{longtable} \r\n  \\end{center}", "\r\n", "\\newpage", "\r\n", sep="")
	cat("\r\n", "\r\n", sep="")
    
    if (i2<immunoassay.environment$n) cat("\n \\newpage \n\n", sep="")
    }
    
	rm(tab,tab1,tab2, i2, i3, pfile)
@

\end{document}