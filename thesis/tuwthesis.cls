\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{tuwthesis}[2026/01/09 TU Wien Thesis Class]

\DeclareOption*{%
  \PassOptionsToClass{\CurrentOption}{scrbook}%
}

\ProcessOptions\relax

% Load scrbook with forwarded options
\LoadClass{scrbook}

%====================================================
% Class Options, Defaults
%====================================================

\newcommand{\supervisor}[1]{\gdef\@supervisor{#1}}
\newcommand{\cosupervisor}[1]{\gdef\@cosupervisor{#1}}
\newcommand{\institute}[1]{\gdef\@institute{#1}}
\newcommand{\university}[1]{\gdef\@university{#1}}
\newcommand{\type}[1]{\gdef\@type{#1}}
\newcommand{\address}[1]{\gdef\@address{#1}}
\newcommand{\acProgram}[1]{\gdef\@acProgram{#1}}
\newcommand{\tulogo}[1]{\gdef\@tulogo{#1}}
\newcommand{\institutelogo}[1]{\gdef\@institutelogo{#1}}

% Default values
\def\@supervisor{\@latex@warning@no@line{No \noexpand\supervisor was given}}
\def\@cosupervisor{\@latex@warning@no@line{No \noexpand\cosupervisor was given}}
\def\@defaultcosupervisor{\@latex@warning@no@line{No \noexpand\cosupervisor was given}}
\def\@institute{\@latex@warning@no@line{No \noexpand\institute was given}}
\def\@university{TU Wien}
\def\@type{BACHELOR THESIS}
\def\@address{}
\def\@acProgram{\@latex@warning@no@line{No \noexpand\acProgram was given}}
\def\@tulogo{preamble/logos/TU_Wien_Logo.pdf}
\def\@institutelogo{preamble/logos/Itp-logo.png}

%====================================================
% Title Page
%====================================================
\renewcommand\maketitle{
\begin{titlepage}
    \pagenumbering{gobble}
    \centering

    % Logos
    \includegraphics[height=72pt]{\@tulogo}
    \hspace{0.25em}
    \includegraphics[height=72pt]{\@institutelogo}
    \\
    \vspace{3em}

    % Type and Title
    \large{\@type} \\
    \vspace{0.5em}
    \Huge{\@title}  \\
    \vspace{2em}

    % Institute & University
    \large{at the \@institute \\ of the \@university} \\
    \vspace{2em}

    % Supervisor
    \large{under supervision of}\\
    \Large{\textbf{\usekomafont{author}\@supervisor}} \\

    % Co-supervisor (optional)
    \ifx\@cosupervisor\@defaultcosupervisor
    \else
        \vspace{2em}
        \large{and}\\
        \Large{\textbf{\usekomafont{author}\@cosupervisor}}\\
    \fi
    \vspace{2em}

    % Author
    \large{by} \\
    \Large{\textbf{\usekomafont{author}\@author}} \\
    \vspace{2em}

    % Address
    \large{\@address}
    \vspace{0.5em}

    % Date
    \large{\@date} \\
    \vspace{4em}

    % Signatures
    \begin{tabular}{p{.47\textwidth}p{.47\textwidth}}
        \rule{6cm}{0.4pt}   & \rule{6cm}{0.4pt} \\
        Signature           & Signature         \\
        Supervisor          & Student           \\
    \end{tabular}

\newpage \ \newpage
\end{titlepage}
}

\linespread{1} % Set linespace

\renewcommand{\frontmatter}{
    \pagenumbering{roman}
    \clearpairofpagestyles
    \ofoot[\pagemark]{\pagemark}
}
\renewcommand{\mainmatter}{
    \pagenumbering{arabic}
    \linespread{1.5} % Set linespace
    \clearpairofpagestyles
    \ihead{\leftmark}
    \ohead{\rightmark}
    \ofoot[\pagemark]{\pagemark}
}

\renewcommand{\backmatter}{
    \cleardoublepage            % start on a new right-hand page
    %\pagenumbering{arabic}      % continue Arabic numbering from main matter
}