\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\ProvidesClass{myposter}[2008/09/18 poster presentation template]
\LoadClassWithOptions{a0poster}

\RequirePackage[left=1cm,top=0.5cm,right=1cm,nohead,nofoot]{geometry}
\RequirePackage{pstricks,pst-node,epsfig}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of some distances
\setlength{\columnsep}{3cm}
\setlength{\columnseprule}{2mm}
\setlength{\parindent}{0.0cm}
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of some variables
\newcommand*{\@inst}{}
\newcommand{\institute}[1]{\gdef\@inst{#1}}
\newcommand*{\@Addr}{}
\newcommand{\address}[1]{\gdef\@Addr{#1}}
\newcommand*{\@Email}{}
\newcommand{\email}[1]{\gdef\@Email{#1}}
\newcommand*{\@Url}{}
\newcommand{\url}[1]{\gdef\@Url{#1}}
\newcommand*{\@Leftlogo}{}
\newcommand{\leftlogo}[1]{\gdef\@Leftlogo{#1}}
\newcommand*{\@Rightlogo}{}
\newcommand{\rightlogo}[1]{\gdef\@Rightlogo{#1}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Background                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\background}{
  \psframe[%
	fillstyle=solid,%
	fillcolor=pksorange%
	%fillstyle=gradient,%
	%gradend=pksblue3,%
	%gradbegin=pksorange,%
	%gradmidpoint=1.
	](0.0\textwidth,0.0\textheight)(1.\textwidth,-1.\textheight)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                Poster                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newenvironment{poster}{
  \vspace*{-0.02\textheight}
  \titlebar
  \begin{minipage}{0.973\textwidth}
  \begin{center}
  \begin{minipage}{0.973\textwidth}
}{
  \end{minipage}
  \end{center}
  \end{minipage}
} 
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                pcolumn                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newenvironment{pcolumn}[1]{
  \begin{minipage}[t][0.85\textheight][b]{#1\textwidth}
  \begin{center}
}{
  \end{center}
  \end{minipage}
}

\newenvironment{pminipage}[1]{
  \hfill
  \begin{minipage}{#1\textwidth}
  %\vspace{1cm}
}{
  \end{minipage}
  \hfill
}

\newenvironment{pminipagef}[2]{
  \hfill
  \begin{minipage}[t][#1]{#2\textwidth}
  \vspace{1cm}
}{
  \end{minipage}
  \hfill
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                pbox                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newrgbcolor{lightblue}{0. 0. 0.80}
\newrgbcolor{white}{1. 1. 1.}
\newrgbcolor{whiteblue}{.80 .80 1.}
\newrgbcolor{pksorange}{0.96 0.54 0.12}
\newrgbcolor{pksblue4}{0 0.2 0.4}
\newrgbcolor{pksblue3}{0.25 0.40 0.55}
\newrgbcolor{pksblue2}{0.65 0.74 0.84}
\newrgbcolor{pksblue1}{0.92 0.94 0.95}

\newcommand{\titlebar}{%
  \vspace*{3.0cm}
  \begin{center}
  \hspace*{-0.01\textwidth}
  \begin{minipage}{0.9\textwidth}
  \psshadowbox[%
	shadowsize=3pt,%
	linewidth=2mm,%
	%framearc=0.1,%
	linearc=1cm,%
	cornersize=absolute,%
	linecolor=pksblue3,%
	fillstyle=solid,%
	fillcolor=white,%
	%fillstyle=gradient,%
	%gradangle=45,%
	%gradbegin=pksorange,%
	%gradend=white,%
	%gradmidpoint=1.0,%
	framesep=1em
	]
	{%
	  \begin{minipage}[t]{0.98\textwidth}
		%%% Logo
		\begin{minipage}[c][8cm][c]{0.1\textwidth}
		  \ifx \@Leftlogo \@empty
		     \hfill
		  \else
			\begin{center}
			\includegraphics[width=0.8\textwidth,angle=0]{\@Leftlogo}
			\end{center}
		  \fi
		\end{minipage}
		%%% Title
		\begin{minipage}[c][8cm][c]{0.79\textwidth}
		  \begin{center}
		  {\sc \huge \@title}\\[10mm]
		  {\large \@author \\[10mm]
		  \ifx \@Email \@empty
			\@inst\\
			\@Addr
		  \else
			\begin{minipage}{0.49\textwidth}
			\centering
			\@inst\\
			\@Addr
			\end{minipage}
			\begin{minipage}{0.49\textwidth}
			\centering
			{\tt \@Email}\\
			\ifx \@Url \@empty \else {\tt \@Url} \fi
			\end{minipage}
		  \fi
		  }
		  \end{center}
		\end{minipage}
		%%% Logo
		\begin{minipage}[c][8cm][c]{0.1\textwidth}
		  \ifx \@Rightlogo \@empty
		  	\hfill
		  \else
			\begin{center}
			\includegraphics[width=0.8\textwidth,angle=0]{\@Rightlogo}
			\end{center}
		  \fi
		\end{minipage}

	  \end{minipage}
	}
  \end{minipage}
  \end{center}
  \vspace*{1cm}
}

\newcommand{\pbox}[1]{%
  \psshadowbox[%
	shadowsize=3pt,%
	linewidth=2mm,%
	%framearc=0.1,%
	linearc=1cm,%
	cornersize=absolute,%
	linecolor=pksblue3,%
	fillstyle=solid,%
	fillcolor=white,%
	%fillstyle=gradient,%
	%gradangle=45,%
	%gradbegin=pksorange,%
	%gradend=white,%
	%gradmidpoint=1.0,%
	framesep=1em]
	{%
	  \begin{minipage}[t]{1.0\textwidth}
#1
	  \end{minipage}
	}\\[1cm]
}

\newcommand{\phighlightbox}[1]{%
  \psshadowbox[%
	shadowsize=3pt,%
	linewidth=2mm,%
	%framearc=0.1,%
	linearc=1cm,%
	cornersize=absolute,%
	linecolor=pksblue4,%
	fillstyle=solid,%
	fillcolor=white,%
	%fillstyle=gradient,%
	%gradangle=45,%
	%gradbegin=pksorange,%
	%gradend=white,%
	%gradmidpoint=1.0,%
	framesep=1em]
	{%
	  \begin{minipage}[t]{1.0\textwidth}
#1
	  \end{minipage}
	}\\[1cm]
}
 
\newcommand{\pboxf}[2]{%
  \psshadowbox[%
	shadowsize=3pt,%
	linewidth=2mm,%
	%framearc=0.1,%
	linearc=1cm,%
	cornersize=absolute,%
	linecolor=pksblue3,%
	fillstyle=solid,%
	fillcolor=white,%
	%fillstyle=gradient,%
	%gradangle=45,%
	%gradbegin=pksorange,%
	%gradend=white,%
	%gradmidpoint=1.0,%
	framesep=1em]{%
	  \begin{minipage}[t][#1]{1.0\textwidth}
#2
	  \end{minipage}
	}\\[1cm]
}
 
\newcommand{\tbox}[2]{%
  \pbox{
  \vspace{-15mm}	
  \psframebox[%
	linewidth=2mm,%
	framearc=0.1,%
	linecolor=pksblue3,%
	fillstyle=solid,%
	fillcolor=pksblue3,%
	framesep=10pt]
	{%
	  \begin{minipage}[t]{0.6\textwidth}
	  \centering
	  \Large \bf
	  \textcolor{white}{#1}
	  \end{minipage}
	} \\[10pt]
#2
  }
}
 
\newcommand{\thighlightbox}[2]{%
  \phighlightbox{
  \vspace{-15mm}	
  \psframebox[%
	linewidth=2mm,%
	framearc=0.1,%
	linecolor=pksblue4,%
	fillstyle=solid,%
	fillcolor=pksblue4,%
	framesep=10pt]
	{%
	  \begin{minipage}[t]{0.6\textwidth}
	  \centering
	  \Large \bf
	  \textcolor{white}{#1}
	  \end{minipage}
	} \\[10pt]
#2
  }
}

\newcommand{\tboxf}[3]{%
  \pboxf{#1}{
  \vspace{-15mm}	
  \psframebox[%
	linewidth=2mm,%
	framearc=0.1,%
	linecolor=pksblue3,%
	fillstyle=solid,%
	fillcolor=pksblue3,%
	framesep=10pt]
	{%
	  \begin{minipage}[t]{0.6\textwidth}
	  \centering
	  \Large \bf
	  \textcolor{white}{#2}
	  \end{minipage}
	} \\[10pt]
#3
  }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                highlight box                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\highlight}[1]{
    \fcolorbox{black}{pksblue1}{
		\begin{minipage}[t]{0.99\textwidth}
		#1
		\end{minipage}
	}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                myfig                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \myfig - replacement for \figure
% necessary, since in multicol-environment 
% \figure won't work

\newcommand{\myfig}[3][0]{
\begin{center}
  \vspace{1.5cm}
  \includegraphics[width=#3\hsize,angle=#1]{#2}
  \nobreak\medskip
\end{center}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                mycaption                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \mycaption - replacement for \caption
% necessary, since in multicol-environment \figure and
% therefore \caption won't work

%\newcounter{figure}
\setcounter{figure}{1}
\newcommand{\mycaption}[1]{
  \vspace{0.5cm}
  \begin{quote}
    {{\sc Fig.} \arabic{figure}: #1}
  \end{quote}
  \vspace{1cm}
  \stepcounter{figure}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                enumerate                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% redefinition of the enumerate counter layout

\renewcommand{\labelenumi}{
  \psset{unit=14pt}
  \begin{pspicture}(0,0)(1,1)
  \psset{fillstyle=solid,fillcolor=black,fillcolor=pksblue3,linecolor=black}
  \rput(0.5,0.5){\pscircle{1}}
  \rput(0.5,0.5){\textcolor{white}{\arabic{enumi}}}
  \end{pspicture}
}
 
\renewcommand{\labelenumii}{
  \textcolor{pksblue3}{\bf \arabic{enumi}.\alph{enumii})}
}
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                itemize                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% redefinition of the itemize counter layout

\renewcommand{\labelitemi}{
  \psset{unit=14pt}
  \begin{pspicture}(0,0)(1,1)
  \psset{fillstyle=solid,fillcolor=black,fillcolor=pksblue3,linecolor=black}
  \rput(0.5,0.5){\pscircle{0.4}}
  \end{pspicture}
}

\renewcommand{\labelitemii}{
  \psset{unit=14pt}
  \begin{pspicture}(0,0)(1,1)
  \psset{fillstyle=solid,fillcolor=black,fillcolor=pksblue3,linecolor=black}
  \rput(0.5,0.5){\pscircle{0.3}}
  \end{pspicture}
}
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                special character             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\RightarrowFilled}{
  \psset{unit=14pt}
  \begin{pspicture}(0,0)(1,1)
  \psset{fillstyle=solid,fillcolor=black,fillcolor=pksblue3,linecolor=black}
  \rput(0.5,0.5){\pscircle{1}}
  \rput(0.5,0.5){\textcolor{white}{\ensuremath{\Rightarrow}}}
  \end{pspicture}
}
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                references                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% redefinition of the references headline

\renewcommand{\refname}{}
