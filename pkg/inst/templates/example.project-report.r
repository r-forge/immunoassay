\documentclass{combine}

\usepackage{combine} 
\usepackage{geometry} 
\usepackage{graphicx}
\usepackage{longtable}
\usepackage{pdflscape}
\usepackage{fancyhdr}
\usepackage[iso]{isodateo}
\usepackage{Sweave}

\geometry{tmargin=2cm,bmargin=2.5cm,lmargin=2cm,rmargin=2cm}
\pagestyle{fancy}
\headheight 25pt

\newcommand{\HRule}{\rule{\linewidth}{0.5mm}}

\begin{document}

\begin{titlepage}

\begin{center}

% Upper part of the page
   
\HRule \\[0.7cm]
\textsc{\huge EXAMPLE Biomarker Study} \\[1.5cm]

\textsc{\huge $A\beta_{1-42}$, $t-tau$ and $p-tau_{181P}$} \\[0.5cm]

\HRule \\[3.5cm]

{ \LARGE \bfseries  Part 1. Project QC summary\\[0.3cm]
                    Part 2. Individual plate reports}\\[1.0cm]

\Large Created in ``R'' and \LaTeX{}

\vfill

% Date
{\large \today}\\[1.0cm]

% Authors - bottom of the page
\begin{minipage}{0.4\textwidth}
\begin{flushleft} \large
\emph{Created by:}\\
Michal J. \textsc{Figurski}, PhD.
\end{flushleft}
\end{minipage}
\begin{minipage}{0.4\textwidth}
\begin{flushright} \large
\emph{Reviewed by:} \\
Michal J. \textsc{Figurski}, PhD.
\end{flushright}
\end{minipage}

\end{center}

\end{titlepage}



\newpage

\mbox{}
\vfill
\begin{center}
\section*{\huge Part 1. Project QC summary}
\end{center}
\vfill

\newpage

\lhead{General QC Summary}
\rhead{Project: \textsc{\Sexpr{immunoassay.environment$project.options$template} }}

<<results=tex, echo=FALSE>>=
    .a  = immunoassay.environment$results
    .a$Type[.a$Type=="Retest"] = "Sample"

    .b  = project.means(.a)
    .b$Kit = factor(.b$Kit)
    .b  = .b[order(.b$Kit, .b$ID, .b$Type, .b$SPL),]
    .bb = list(c("St7"),c("St6","St5"),c("St4","St3"),c("St2","St1"),c("pool48","pool52"))

    # Project-level plots:
    for (i2 in 1:length(.bb)) {
        CairoPDF(file=paste(immunoassay.environment$path, "/Tex/pdf/psum-", i2, sep=""), width=14, height=8)
            par(mfrow=c(1,2))
            .bc = c(1, length(unique(.b$ID)))
            if (i2 == 1)     .bc = c(1,15,17,20,24,35,length(unique(.b$ID)))
            if (i2 %in% 2:4) .bc = c(1,17,20,24,35,length(unique(.b$ID)))
            project.qcplot(.b, samples=.bb[[i2]], breaks=.bc)
        dev.off()
        
        cat("\\vfill\n")
        cat("\\begin{center}\n", sep="")
        cat("   \\includegraphics[width=1.01\\textwidth]{\"", "pdf/psum-",i2, "\"}\n", sep="")
        cat("\\end{center}\n\n", sep="")
        cat("\\vfill\n")
        }
@

\newpage

\lhead{ }
\rhead{ }

\mbox{}
\vfill
\begin{center}
\section*{\huge Part 2. Individual plate reports}
\end{center}
\vfill

\newpage

\lhead{Data interpretation report}
\rhead{Acquisition Time: \Sexpr{attr(immunoassay.environment$l.run,"Date")} }
\rfoot{Page: \thepage}
\cfoot{Report Date: \today}

\setcounter{secnumdepth}{-1 }

\pagestyle{combine}
\begin{papers}

<<results=tex, echo=FALSE>>=
	for (i1 in immunoassay.environment$files) {
        i2  = substr(i1,0, nchar(i1)-4)
		cat(paste("\\import{\"",i2,"\"} \n\n", sep=""))
		}
        
    # Clean-up
    rm(.a,.b,.bb,.bc,i1,i2)
@

\end{papers}
\end{document}
