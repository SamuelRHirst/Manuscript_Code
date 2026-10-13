# Collapse RepeatMasker / GraffiTE labels from the Catharus library into broad repeat classes
function cls(x,   y) {
    y = tolower(x)
    if (y == "simple_repeat")                       return "Simple repeat"
    if (y == "low_complexity")                      return "Low complexity"
    if (y ~ /^satellite/)                           return "Satellite"
    if (y ~ /^(rrna|trna|snrna|scrna|srprna)$/)     return "Small RNA"
    if (y ~ /^sine/)                                return "SINE"
    if (y ~ /^line/ || y ~ /^ple(\/|$)/)            return "LINE"
    if (y ~ /^ltr/)                                 return "LTR"
    if (y ~ /helitron/ || y ~ /^rc(\/|$)/)          return "Rolling circle"
    if (y ~ /^dna(\/|$)/ || y ~ /^mite(\/|$)/)      return "DNA"
    if (y ~ /^(unknown|unspecified)/)               return "Unknown"
    return "Other"
}
