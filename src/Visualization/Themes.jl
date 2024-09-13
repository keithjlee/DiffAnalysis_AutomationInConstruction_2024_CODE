# text sizes
tickFontSize = 12
labelFontSize = 14
titleFontSize = 14

# font
typeface = "IBM Plex Mono" #IBM Plex Sans, Libre Caslon Text, Times New Roman, Neue Haas Unica Pro

aic = Theme(
    backgroundcolor = :white,

    palette = (color = [kjl_blue, kjl_green, kjl_pink, kjl_orange, kjl_gray, kjl_darkgray, :black],
        marker = [:circle, :rect, :utriangle, :xcross]),

    Axis = (aspect = 1,
        backgroundcolor = :transparent,
        xlabelfont = typeface,
        xlabelsize = labelFontSize,
        ylabelfont = typeface,
        ylabelsize = labelFontSize,
        xticklabelfont= typeface,
        xticklabelsize = tickFontSize,
        yticklabelfont = typeface,
        yticklabelsize = tickFontSize,
        titlefont = typeface,
        titlesize = titleFontSize,
        xgridcolor = kjl_gray,
        ygridcolor = kjl_gray,
        xgridvisible = false,
        ygridvisible = false,
        xtickalign = 1,
        ytickalign = 1,
        xticksmirrored = true,
        yticksmirrored = true
        ),

    Axis3 = (backgroundcolor = :transparent,
        aspect = (1,1,1),
        xgridcolor = kjl_gray,
        ygridcolor = kjl_gray,
        zgridcolor = kjl_gray,
        titlefont = typeface,
        xlabelfont = typeface,
        ylabelfont = typeface,
        zlabelfont = typeface,
        xticklabelfont = typeface,
        yticklabelfont = typeface,
        zticklabelfont = typeface,
        titlesize = titleFontSize,
        xlabelsize = labelFontSize,
        ylabelsize = labelFontSize,
        zlabelsize = labelFontSize,
        xticklabelsize = tickFontSize,
        yticklabelsize = tickFontSize,
        zticklabelsize = tickFontSize,
        azimuth = -3pi / 4,
        elevation = pi/8
        ),

    Lines = (linewidth = 2,),

    LineSegments = (linewidth = 2,),

    Colorbar = (labelfont = typeface,
        labelsize = labelFontSize,
        ticksvisible = false,
        spinewidth = 0,
        ticklabelfont = typeface,
        ticklabelsize = tickFontSize,
        colormap = trans2black,),

    Legend = (backgroundcolor = :white,
        framecolor = :black,
        labelsize = tickFontSize,
        labelfont = typeface,
        titlefont = typeface),

    Heatmap = (colormap = trans2black,),

    Surface = (colormap = trans2black,),

    Spy = (colormap = trans2black,),

    Scatter = (strokewidth = 1,
        cycle = [:color, :marker],
        strokecolor = :black),

    Hist = (color = kjl_blue,
        strokecolor = :black),

    BarPlot = (color = kjl_blue,
        gap = 0,
        strokewidth = 0,
        strokecolor = :white,),    

    Text = (font = typeface,),
)