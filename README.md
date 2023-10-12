# HIL-FaultMitigation

Welcome to the **HIL-FaultMitigation** repository: Code to perform Hardware-In-Loop (HIL) activities on engines using INCA and PUMA data acquisition systems. 

We use the HIL system to mitigate pressure shift faults on a 7.6-liter Navistar diesel engine. The code primarily involves reading engine data from INCA and PUMA, preprocessing it, implementing the controller to maintain torque profiles or mitigate faults, and finally, writing the engine inputs to the ECU.

\begin{table}[]
\setlength{\tabcolsep}{2pt}
\renewcommand{\arraystretch}{1.2}
\fontsize{9pt}{9pt}\selectfont
\centering
\caption{Hardware and API details for INCA and PUMA}
\label{tab:HIL_details}
\resizebox{255}{60}{%
\begin{tabular}{|l|l|l|}
\hline
\rowcolor[HTML]{BBE0E3} 
  & \multicolumn{1}{c|}{\cellcolor[HTML]{BBE0E3}\textbf{INCA}}                 & \multicolumn{1}{c|}{\cellcolor[HTML]{BBE0E3}\textbf{PUMA}}           \\ \hline
\rowcolor[HTML]{F3F9FA} 
  \textbf{Function} &
  \begin{tabular}[c]{@{}l@{}}Read ECU data and \\ write the control signals\end{tabular} &
  \begin{tabular}[c]{@{}l@{}}Read time-based sensor data \\ (such as mass flow, temperature, \\  pressure, torque, and speed)\end{tabular} \\ \hline
\rowcolor[HTML]{E7F3F4} 
 \textbf{GUI   Name}       & INCA V7.2.17                  & Kvaser CANKing          \\  \hline 
\rowcolor[HTML]{F3F9FA} 
 \textbf{Adapter   Name}   & ES581.3                       & Kvaser Leaf Light V2    \\ \hline
\rowcolor[HTML]{E7F3F4} 
 \textbf{Adapter   Image}  & \multicolumn{1}{c|}{\cellcolor[HTML]{E7F3F4}{\includegraphics[width=0.12\textwidth, height=8.5mm]{Figures/ETAS_Adapter.png} } } &   \multicolumn{1}{c|}{\cellcolor[HTML]{E7F3F4}{\includegraphics[width=0.12\textwidth, height=8.5mm]{Figures/Kvaser_Adapter.png} } }            \\ \hline
\rowcolor[HTML]{F3F9FA} 
 \textbf{Python   Library} & INCA-python                   & CANlib and cantools     \\ \hline
\rowcolor[HTML]{E7F3F4} 
\textbf{Sampling   Rate}  & \begin{tabular}[c]{@{}l@{}} Synchronous with \\ engine speed\end{tabular} & ½, 1, 2, 4, 5, 10 Hertz \\ \hline
\end{tabular}%
}
\end{table}
