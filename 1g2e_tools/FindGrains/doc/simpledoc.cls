%==============================
% simpledoc.cls
% miniclasse to write very simple docs  
% Vincent Richefeu 2003, 2006
% last modifications : 15/02/2006
%==============================
\ProvidesClass{gdmdoc}
\LoadClass[twoside,a4paper]{refart}

%\RequirePackage[french]{varioref}
%\RequirePackage[french]{babel} 
%\RequirePackage[latin1]{inputenc}  % for PC
\RequirePackage[applemac]{inputenc} % for apple

\RequirePackage{amsmath}
\RequirePackage{graphicx}
\RequirePackage{float}
\RequirePackage[sort]{natbib}
\RequirePackage{tocloft}

%
\RequirePackage{listings}
\RequirePackage{color}


\definecolor{myblue}{rgb}{0.0,0.53333,0.80784}% RGB #0088CE, Pantone Process Blue
\definecolor{myblue2}{rgb}{0.0,0.2,0.6}
\definecolor{myyellow}{rgb}{1.0,1.0,0.95}

\lstdefinestyle{cppCode}{
language = C++,
rulecolor=\color{myblue2},
backgroundcolor=\color{myyellow},
frame = single,
numbers = left,
breaklines = true,
numberstyle = \small\ttfamily\color{myblue},
basicstyle = \small\ttfamily\color{myblue},
texcl = true,
commentstyle = \color[rgb]{0.1,0.1,0.1},
keywordstyle = \small\ttfamily\color{myblue2},
}

\lstdefinestyle{simfile}{
language = C++,
rulecolor=\color{myblue2},
backgroundcolor=\color{myyellow},
frame = single,
numberstyle=\small\ttfamily\color{myblue},
basicstyle=\small\ttfamily\color{myblue},
breaklines=true,
emph = {Comment,Simulation,Physic,Data_table,Driver,Properties,Spy,Informer,Sample,Network} , emphstyle=\bfseries\color{myblue2}
}

%%
