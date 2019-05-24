%% file ethbrief3.cls
%% Version 3.1 new form of corporated design -- © Dieter Hennig
%% based an Version 2 and with the helpfull sugesstions
%% for Version 3.1 by Niklas Beisert

%% --- Identification part ---
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{ethbrief3}[2013/07/04 ETH Brief Version 3.1]

%% --- Package loading ---
\RequirePackage{graphics}
\RequirePackage{ifthen}

%% --- use of hyperref only if present --- 
\newcommand{\eth@href}[2]{\@ifundefined{Hy@Warning}{#2}{\href{#1}{#2}}}

% Multiple language support (Robert Mach, Dieter Hennig)
\newboolean{@ngerman} \setboolean{@ngerman}{false}
\newboolean{@english} \setboolean{@english}{false}
\newboolean{@EthUzhLogo} \setboolean{@EthUzhLogo}{false}
\newboolean{@DLogo} \setboolean{@DLogo}{false}

%% --- Declaration of options ---
% Multiple language support (Robert Mach et al.)
\DeclareOption{german}{\setboolean{@ngerman}{true}}
\DeclareOption{ngerman}{\setboolean{@ngerman}{true}}
\DeclareOption{english}{\setboolean{@english}{true}}
\DeclareOption{ETH}{\setboolean{@EthUzhLogo}{false}}    % ETH Logo (default)
\DeclareOption{ETHUZH}{\setboolean{@EthUzhLogo}{true}}  % ETH-UZH Logo
\DeclareOption{DepLogo}{\setboolean{@DLogo}{true}}      % Departement-Logo
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{letter}}


%% --- Execution of options ---
\ProcessOptions
\LoadClass[a4paper]{letter}


%% --- Special commands for package and class files ---
\usepackage[T1]{fontenc}
\usepackage[pagecontinue=true]{pageslts}
\pagenumbering{arabic}
%\usepackage{microtype}

%%%%%%%%%%%%%%%%%%
% With Helvetica fonts
%%%%%%%%%%%%%%%%%%

% Use the following lines if the hevetica fonts should be used
% (Frank Schweitzer, Robert Mach, Dieter Hennig, Niklas Beisert)
\newfont{\eth@xiiiv}{phvr8t at 8.5pt}     % header-font regular
\newfont{\eth@xiiivb}{phvb8t at 8.6pt}    % header-font bold

%\renewcommand{\familydefault}{phv} % Helvetica 1.methode
%\renewcommand{\sfdefault}{phv}     % Helvetica 1.methode
\usepackage[scaled=.95]{helvet}    % Helvetica 2.methode
%\usepackage[scaled=.95]{uarial}    %  URW A030
%\renewcommand{\familydefault}{uhv} % Helvetica UWR Nimbus Sans L 3.methode
%\renewcommand{\sfdefault}{uhv}     % Helvetica UWR Nimbus Sans L 3.methode

\renewcommand{\familydefault}{\sfdefault} % all methods

%% --- Text dimensions  ---
\setlength{\textwidth}{166mm}
\setlength{\textheight}{237mm}
\setlength{\voffset}{-4.9mm}
\setlength{\topmargin}{0mm}
\setlength{\headheight}{0mm}
\setlength{\headsep}{5mm}
\setlength{\oddsidemargin}{-3mm}
\setlength{\evensidemargin}{\oddsidemargin}
\setlength{\hoffset}{0mm}
\setlength{\topskip}{30mm}
\setlength{\footskip}{24mm}

%% Enclosed, Carbon Copy, ...
\renewcommand{\encl}[1]{%
\par\noindent
\parbox[t]{\textwidth}{
\@hangfrom{\ethenclname }
\ignorespaces #1\strut}}%

%% --- new declaration of language renaming ---
\renewcommand{\cc}{\ethccname\ }
\renewcommand{\ps}{\ethpsname\ }

%% --- include of the signature file ---
\newcommand{\eth@signature}{}
\newcommand{\eth@signaturefile}{}
\newcommand{\signaturefile}[2][0.6cm]{\renewcommand{\eth@signaturefile}{\raisebox{#1}{#2}}}
\renewcommand{\signature}[1]{\renewcommand{\eth@signature}{\makebox[0cm][l]{\smash{\eth@signaturefile}}#1}}
%% --- left alignment of text, interlinepenalty
\raggedright
\interlinepenalty=200
%% --- work around \begin{letter} for serial letters (Niklas Beisert)
\newcommand{\recipient}[1]{\@processto{\leavevmode\ignorespaces #1}}

% Multiple language support (Robert Mach, Niklas Beisert)
\newcommand{\EthLang}[2]{\ifthenelse{\boolean{@ngerman}}{#1}{#2}}

%% --- simplified language support ---
%% --- internal CD-names, not switched by babel
\newcommand{\ethpagename}{\EthLang{Seite}{page}}
\newcommand{\ethenclname}{\EthLang{Beilage(n)}{Enclosure(s)}}
\newcommand{\ethccname}{\EthLang{Kopie}{CC}}
\newcommand{\ethpsname}{\EthLang{PS}{PS}}
\newcommand{\ethphonename}{\EthLang{Telefon}{Phone}}
\newcommand{\ethfaxname}{\EthLang{Telefax}{Fax}}
\newcommand{\ethmobilename}{\EthLang{Mobil}{Mobile}}
%% universal location names
\newcommand{\Zurich}{\EthLang{Z\"urich}{Zurich}}
\newcommand{\Switzerland}{\EthLang{Schweiz}{Switzerland}}

%% --- Header definitions for first page (Robert Mach, Dieter Hennig) ---
\renewcommand{\ps@firstpage}
{
 \renewcommand{\@oddhead}{}%
 \ifthenelse{\boolean{@DLogo}} % Departement/Institut-Logo for the footer
    {% work with Dep-Logo on first page for ETH only page
     % version 3 foot has space for departement-logo or is empty
     \ifthenelse{\boolean{@EthUzhLogo}}{ % with ETH-UZH Logo no third logo 
                 \renewcommand{\@oddfoot}{}}%
                {% ETH Logo
                 \renewcommand{\@oddfoot}{\mbox{\DLogo[90mm]}}}
      }% end ETH-UZH
    {% without Departement/Institut-Logo
      \renewcommand{\@oddfoot}{}
 }% finish Departement/Institut-Logo for the footer
 \renewcommand{\@evenhead}{}\renewcommand{\@evenfoot}{}
} % end firstpage
%% --- and following pages --
\renewcommand{\ps@headings}
{
 \ifthenelse{\boolean{@EthUzhLogo}}
 {% ETH-UZH Logo
    \renewcommand{\@oddhead}{\raisebox{-17mm}{\ifx\@empty\eth@SubjectTextHead{\ } \else\textbf{\eth@xiiivb\eth@SubjectTextHead}\fi}}%
 }
 {% ETH Logo
    \renewcommand{\@oddhead}{\mbox{\ifx\@empty\eth@SubjectTextHead{\ } \else\textbf{\eth@xiiivb\eth@SubjectTextHead}\fi}}%
 }
 \renewcommand{\@oddfoot}%
 {% uncomment the next line, in case of trouble with page-numbering for serial letters
 \hfill\raggedleft\mbox{\eth@xiiiv\ethpagename~\thepage~/~\lastpageref*{pagesLTS.arabic}}
 }
 \renewcommand{\@evenhead}{\@oddhead}
 \renewcommand{\@evenfoot}{\@oddfoot}
} % end all other pages

\pagestyle{headings}

%% --- internal text parameters (Niklas Beisert) 
\newcommand{\eth@FirstUnitNameText}{}
\newcommand{\EthFirstUnitName}[1]{\renewcommand{\eth@FirstUnitNameText}{#1}}
\newcommand{\eth@SecondUnitNameText}{}
\newcommand{\EthSecondUnitName}[1]{\renewcommand{\eth@SecondUnitNameText}{\\#1}}
\newcommand{\eth@PersNameText}{}
\newcommand{\EthPersName}[1]{\renewcommand{\eth@PersNameText}{\\#1}}
\newcommand{\eth@PersPositionText}{}
\newcommand{\EthPersPosition}[1]{\renewcommand{\eth@PersPositionText}{\\#1}}
\newcommand{\eth@NameText}{}
\newcommand{\EthName}[1]{\renewcommand{\eth@NameText}{#1}}
\newcommand{\eth@AddressText}{\eth@OfficeText\eth@StreetText\eth@TownText\eth@CountryText}
\newcommand{\EthAddress}[1]{\renewcommand{\eth@AddressText}{\\#1}}
\newcommand{\eth@PersInfoText}{}
\newcommand{\EthPersInfo}[1]{\renewcommand{\eth@PersInfoText}{\\#1}}
\newcommand{\eth@OfficeText}{}
\newcommand{\EthOffice}[1]{\renewcommand{\eth@OfficeText}{\\#1}}
\newcommand{\eth@StreetText}{}
\newcommand{\EthStreet}[1]{\renewcommand{\eth@StreetText}{\\#1}}
\newcommand{\eth@TownText}{}
\newcommand{\EthTown}[1]{\renewcommand{\eth@TownText}{\\#1}}
\newcommand{\eth@CountryText}{}
\newcommand{\EthCountry}[1]{\renewcommand{\eth@CountryText}{, #1}}
%% --- can add several phone numbers
\newcommand{\eth@PhoneText}{}
\newcommand{\EthPhone}[1]{\g@addto@macro\eth@PhoneText{\\\makebox[3.5em][l]{\ethphonename}#1}}
\newcommand{\eth@FaxText}{}
\newcommand{\EthFax}[1]{\g@addto@macro\eth@FaxText{\\\makebox[3.5em][l]{\ethfaxname}#1}}
\newcommand{\eth@MobileText}{}
\newcommand{\EthMobile}[1]{\g@addto@macro\eth@MobileText{\\\makebox[3.5em][l]{\ethmobilename}#1}}
\newcommand{\eth@EmailText}{}
\newcommand{\EthEmail}[1]{\g@addto@macro\eth@EmailText{\\\eth@href{mailto:#1}{#1}}}
\newcommand{\eth@WebText}{}
\newcommand{\EthWeb}[1]{\g@addto@macro\eth@WebText{\\\eth@href{http://#1}{#1}}}

\newcommand{\eth@SubjectText}{\vspace{-6mm}}
\newcommand{\eth@SubjectTextHead}{\vspace{-6mm}}
\newcommand{\subject}[2][]{\renewcommand{\eth@SubjectText}{#2}%
  \ifthenelse{\equal{#1}{}}%
  {\renewcommand{\eth@SubjectTextHead}{#2}}%
  {\renewcommand{\eth@SubjectTextHead}{#1}}%
 }

\newcommand{\eth@LocationText}{\Zurich}
\renewcommand{\location}[1]{\renewcommand{\eth@LocationText}{#1}}

%% --- backward compability to v3.0
\newcommand{\EthAdress}{\EthAddress}
\newcommand{\Betreff}{\subject}
\newcommand{\EthLocation}{\location}

%% --- compability to the letter class
\renewcommand{\address}{\EthAddress}
\renewcommand{\name}{\EthPersName}
\renewcommand{\telephone}{\EthPhone}

%% --- Replacement for the older ETHlogo.sty (Dieter Hennig)
\ifthenelse{\boolean{@EthUzhLogo}}
  { % if we want a letter with ETH-UZH logo, there are different Logos in German and English for UZH
    \ifthenelse{\boolean{@ngerman}}
      {% german logo
   	     \newcommand{\Logo}[1][\textwidth]{\resizebox{#1}{!}{\includegraphics{ETH-UZHlogo_13}}}
      }% end german logo	
      {% english logo
  	     \newcommand{\Logo}[1][\textwidth]{\resizebox{#1}{!}{\includegraphics{ETH-UZHlogo-en_13}}}
      }% end english logo
  }% end ETH/ETH-UZH logo
  { % else the ETH logo is used
        \newcommand{\Logo}[1][\textwidth]{\resizebox{#1}{!}{\includegraphics{Logo_FAU_DinA4_RGB}}}
  }

% department logo
\newcommand{\eth@DepLogoFile}{Deplogo_13}
\newcommand{\EthDepLogoFile}[1]{\renewcommand{\eth@DepLogoFile}{#1}}
\newcommand{\DLogo}[1][textwidth]{\resizebox{#1}{!}{\includegraphics{\eth@DepLogoFile}}}

% letter head
\renewcommand{\opening}[1]{%
	\thispagestyle{firstpage} % Special style for firstpage
  	\setlength{\unitlength}{1mm}
  %% The picture environment is defined in this way
  %% to ensure that pages following page 1 are right
  %% (Robert Mach, Dieter Hennig)
  	\begin{picture}(166,63)
    	\put(-1,88){
      	% Second picture environment to set origin of
      	% coordinate system to top left.
      	\begin{picture}(0,0)
      	   \ifthenelse{\boolean{@EthUzhLogo}} % Dieter Hennig
          { % letter with ETH-UZH logo
        	    \put(0,-8){\makebox(0,0)[lt]{\Logo[208.3mm]}}
        	    \put(0,-53){\makebox(0,0)[lt]{%
          	\begin{minipage}[t]{98mm}
            	\raggedright \toname \\ \toaddress \par
          	\end{minipage}
       	    }}
            \put(98,-34){\makebox(0,0)[lt]{%
            \begin{minipage}[t]{68mm}
                \eth@xiiiv
                \baselineskip=11pt
                \eth@NameText%
                \eth@PersNameText%
                \eth@PersPositionText\\
                \eth@FirstUnitNameText%
                \eth@SecondUnitNameText%
                \eth@AddressText\\~
                \eth@PhoneText%
                \eth@FaxText%
                \eth@MobileText%
                \eth@EmailText%
                \eth@WebText%
                \eth@PersInfoText%
            \end{minipage}
            }}
          } % end ETH-UZH logo
          { % letter only with the ETH logo
        	       \put(0,-16){\makebox(0,0)[lt]{\Logo[52mm]}}
        	       \put(0,-53){\makebox(0,0)[lt]{%
          	   \begin{minipage}[t]{98mm}
            	   \raggedright \toname \\ \toaddress \par
          	   \end{minipage}
       	       }}
               \put(98,-18){\makebox(0,0)[lt]{%
               \begin{minipage}[t]{68mm}
                   \eth@xiiivb
                   \baselineskip=11pt
                   \eth@FirstUnitNameText%
                   \eth@SecondUnitNameText\\~\\
                   \eth@xiiiv
                   \baselineskip=11pt
                   \eth@NameText%
                   \eth@PersNameText%
                   \eth@PersPositionText%
                   \eth@AddressText\\~
                   \eth@PhoneText%
                   \eth@FaxText%
                   \eth@MobileText%
                   \eth@EmailText%
                   \eth@WebText%
                   \eth@PersInfoText%
              \end{minipage}
            }}
          } % end ETH logo
      \end{picture}
    }
  \end{picture}
\par
\vspace{-24mm}
\eth@LocationText, \@date
\par
\vspace{5mm}
\ifx\@empty\eth@SubjectText\else\textbf{\eth@SubjectText}\vspace{5mm}\fi
\par{#1}
\par
\flushbottom
} % end \renewcommand{\opening}...

% signature
\renewcommand{\closing}[1]{\par \ \par \nobreak\stopbreaks\noindent%
  \parbox{\indentedwidth}{\raggedright\ignorespaces  #1\\[6\medskipamount]
%% --- different fallback-cases for \signature 
  \ifx\@empty\eth@signature\eth@PersNameText\else\eth@signature\fi\strut}%
  \ \par
  \ \par
  }
%% end of file  
