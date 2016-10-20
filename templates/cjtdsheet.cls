% cjtdsheet.cls
%
% Datasheet class used by the CompliantJointToolbox
%
% 2016-Oct-02 Jörn Malzahn
% --- Class structure: identification part
% ---
\ProvidesClass{cjtdsheet}[2016/10/02 version 0.1 CompliantJointToolbox Datasheet]
\NeedsTeXFormat{LaTeX2e}

% --- Initial code part
% --- 
%    Loads packages used internally by the class file
%\InputIfFileExists{cjtdsheet.cfg}
\input{cjtdsheet.cfg}

% --- Class structure: declaration of options part
% ---
% This class extends the article class. Read all the documentclass options; pass them to article,
% unless the file "<currentoption>.cjtdsheet" exists, then it is loaded

\DeclareOption*{\InputIfFileExists{\CurrentOption.cjtdsheet}{}{%
\PassOptionsToClass{\CurrentOption}{article}}}



% --- Class structure: execution of options part
% ---
\ProcessOptions \relax

% --- Class structure: declaration of options part
% ---
\LoadClass{article}

% Page layout
\RequirePackage[left=1.5cm,right=1.5cm,top=1.5cm,bottom=1.5cm]{geometry}

\RequirePackage{multicol}

% Fonts
\RequirePackage[sc]{mathpazo}       % Palatino (not Palatino Linotype!) as font with serifs
\linespread{1.05}                   % Palatino needs more space between lines
\RequirePackage[scaled=.90]{helvet} % Helvetica as font without serifs
\RequirePackage{courier} 			% Courier as type writer font

\RequirePackage{fancyhdr}
\fancypagestyle{cjtpage}{%
\fancyhf{} % clear all six fields
\renewcommand{\headrulewidth}{4pt}
\renewcommand{\footrulewidth}{0pt}
\fancyfoot[L]{Created with the Compliant Joint Toolbox}
\fancyfoot[R]{\today}
}
\pagestyle{cjtpage}



\RequirePackage{graphicx} % for pdf, bitmapped graphics files
\RequirePackage{amsmath} % assumes amsmath package installed
\RequirePackage{amssymb}  % assumes amsmath package installed


\RequirePackage{tabularx}
\RequirePackage{booktabs}
\RequirePackage{color}
\RequirePackage{colortbl}
%% Define user colors using the RGB model
\definecolor{black}{rgb}{0.0,0.0,0.0}
\definecolor{darkdarkgray}{rgb}{0.6,0.6,0.6}
\definecolor{darkgray}{rgb}{0.8,0.8,0.8}
\definecolor{lightgray}{rgb}{0.95,0.95,0.95}
\definecolor{white}{rgb}{1.0,1.0,1.0}
\definecolor{cjtblue}{rgb}{ 0.0314, 0.2157, 0.3882}
\definecolor{cjtlightblue}{rgb}{ 0.4314, 0.6157, 0.7882}
\definecolor{cjtred}{rgb}{ 0.7529, 0, 0}


\RequirePackage[justification=raggedright,
                singlelinecheck=false, 
                font={bf,color=cjtred},
                labelfont+={it},
                position=top,
                skip=0.5\baselineskip
                ]{caption}
\RequirePackage{subcaption}


\newcommand{\actuatorname}[1]
{
\begin{Huge}
	\Huge{\textcolor{cjtred}{\textbf{#1}}}
\end{Huge}

\vspace{0.3cm}
\noindent\rule{\textwidth}{4pt}
}

