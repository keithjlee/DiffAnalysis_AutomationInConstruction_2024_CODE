# base colours
const kjl_pink = colorant"#f52565"
const kjl_blue = colorant"#2559f5"
const kjl_green = colorant"#58bc82"
const kjl_turquoise = colorant"#119da4"
const kjl_orange = colorant"#f5bb25"
const kjl_tan = colorant"#e9d985"
const kjl_brown = colorant"#c3b299"
const kjl_gray = colorant"#cbd4c2"
const kjl_darkgray = colorant"#50514f"

# Legacy
const caitlin_blue = colorant"#3EA8DE"
const caitlin_pink = colorant"#FF7BAC"

# Colour gradients

const pink2blue = cgrad([kjl_pink, :white, kjl_blue])

const white2blue = cgrad([:white, kjl_blue])
const white2pink = cgrad([:white, kjl_pink])
const white2black = cgrad([:white, :black])

const trans2blue = cgrad([:transparent, kjl_blue])
const trans2pink = cgrad([:transparent, kjl_pink])
const trans2black = cgrad([:transparent, :black])
const trans2white = cgrad([:transparent, :white])