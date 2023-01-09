function spec = plotSpec()
%PLOTSPEC Summary of this function goes here
%   Detailed explanation goes here
lineWidth = 1.5;
fontSize  = 11;


anklePlotSpec       = [{'-^'},{'Color'},{'#0072BD'},{'LineWidth'},{lineWidth}];
hipPlotSpec         = [{'-o'},{'Color'},{'#D95319'},{'LineWidth'},{lineWidth}];
toePlotSpec         = [{'-s'},{'Color'},{'#7E2F8E'},{'LineWidth'},{lineWidth}];

LIP_PlotSpec        = [{'Color'},{'#92C1E6'},{'LineWidth'},{lineWidth}];
LIPPFW_PlotSpec     = [{'Color'},{'#FF9F3A'},{'LineWidth'},{lineWidth}];
VHIP_PlotSpec       = [{'Color'},{'#CA78DB'},{'LineWidth'},{lineWidth}];
VHIPPFW_PlotSpec    = [{'Color'},{'#79AD00'},{'LineWidth'},{lineWidth}];

statePlotSpec       =[{'Color'},{'#B1B1B1'},{'LineWidth'},{lineWidth};
                      {'Color'},{'#EDB120'},{'LineWidth'},{lineWidth}];
                  
controlPlotSpec     =[{'Color'},{'#0A0708'},{'LineWidth'},{lineWidth};
                      {'Color'},{'#970C10'},{'LineWidth'},{lineWidth}];
                 
latexFormatting     =[{'interpreter'},{'latex'},{'FontSize'},{fontSize}];                  
                  

spec.LIP        = LIP_PlotSpec;
spec.LIPPFW     = LIPPFW_PlotSpec;
spec.VHIP       = VHIP_PlotSpec;
spec.VHIPPFW    = VHIPPFW_PlotSpec;

spec.ankle      = anklePlotSpec;
spec.hip        = hipPlotSpec;
spec.toe        = toePlotSpec;

spec.state      = statePlotSpec;
spec.control    = controlPlotSpec;

spec.ltxFMT     = latexFormatting;

end

