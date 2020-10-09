---
title: "Data-Preparation-Workflow"
output: 
  html_document:
    fig_width: 10
    fig_height: 8
vignette: >
  %\VignetteIndexEntry{Data-Preparation-Workflow}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
## load the library
library(DEE2HsapienData)

## Creat a 'temp' directory in current working directory if it is not created
## already. In this tutorial, we use it as the workspace 
tempDirPath <- paste(getwd(), 'temp', sep = '/')
dir.create(tempDirPath)
#> Warning in dir.create(tempDirPath): '/home/siyuan/Repos/DEE2HsapienData/
#> vignettes/temp' already exists

## get the 'list' object. 
dee2Data <- GetDEE2Data(
 system.file("extdata", "SelectedDEE2_v2.ods", package = "DEE2HsapienData"),
 paste(tempDirPath, 'DEE2Data', sep = '/') , 
 'hsapiens')
#> Parsed with column specification:
#> cols(
#>   .default = col_character(),
#>   index = col_double(),
#>   ReleaseDate = col_datetime(format = ""),
#>   LoadDate = col_datetime(format = ""),
#>   spots = col_double(),
#>   bases = col_double(),
#>   spots_with_mates = col_double(),
#>   avgLength = col_double(),
#>   size_MB = col_double(),
#>   InsertSize = col_double(),
#>   InsertDev = col_double(),
#>   Study_Pubmed_id = col_double(),
#>   ProjectID = col_double(),
#>   TaxID = col_double(),
#>   g1k_pop_code = col_logical(),
#>   source = col_logical(),
#>   g1k_analysis_group = col_logical(),
#>   Subject_ID = col_logical(),
#>   Disease = col_logical(),
#>   Affection_Status = col_logical(),
#>   Analyte_Type = col_logical()
#>   # ... with 3 more columns
#> )
#> See spec(...) for full column specifications.
#> For more information about DEE2 QC metrics, visit
#>     https://github.com/markziemann/dee2/blob/master/qc/qc_metrics.md

## This operation may requires up to 8GB of memory
dee2Data <- tximportDee2Data('hsapiens', dee2Data, )
#> [1] "Downloading abundance files from DEE2 web API for each SRR accessions,\n        This step could take a very long time if the number of SRR accessions is\n        huge."
#> Note: importing `abundance.h5` is typically faster than `abundance.tsv`
#> reading in files with read_tsv
#> 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364 365 366 367 368 369 370 371 372 373 374 375 376 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 392 393 394 395 396 397 398 399 400 401 402 403 404 405 406 407 408 409 410 411 412 413 414 415 416 417 418 419 420 421 422 423 424 425 426 427 428 429 430 431 432 433 434 435 436 437 438 439 440 441 442 443 444 445 
#> summarizing abundance
#> summarizing counts
#> summarizing length

## Normalize gene-level count data with Deseq2 normalization
dESeq2NormalizedData <- Deseq2Normalization(dee2Data)
## Normalize gene-level count data with Deseq2 normalization
dESeq2NormalizedDataFromTxLevelData <- Deseq2NormalizationFromTranscriptLevelData(dee2Data)
#> using counts and average transcript lengths from tximport
#> using 'avgTxLength' from assays(dds), correcting for library size
## Normalize gene-level count data with Deseq2 normalization
tmmNormalizedData <- TMMNormalization(dee2Data)
#> Warning in filterByExpr.DGEList(dgeList): All samples appear to belong to the
#> same group.
## Normalize gene-level count data with Deseq2 normalization
tmmNormalizedDataFromTxLevelData <- TMMNormalizationFromTranscriptLevelData(dee2Data)
#> Warning in filterByExpr.DGEList(y): All samples appear to belong to the same
#> group.
#> Warning in normarg_mcols(value, class(x), length(x)): You supplied metadata columns of length 13800 to set on an object of
#>   length 14603. However please note that the latter is not a multiple of
#>   the former.

PcaPlot(dESeq2NormalizedData$normalizedSExpr)
#> Warning in if (plot == "2d") {: the condition has length > 1 and only the first
#> element will be used
#> Warning: `arrange_()` is deprecated as of dplyr 0.7.0.
#> Please use `arrange()` instead.
#> See vignette('programming') for more help
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_warnings()` to see where this warning was generated.
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
#> Warning in normalizePath(f2): path[1]="webshot16db861d8feb3.png": No such file
#> or directory
#> Warning in file(con, "rb"): cannot open file 'webshot16db861d8feb3.png': No such
#> file or directory
#> Error in file(con, "rb"): cannot open the connection
PcaPlot(dESeq2NormalizedData$normalizedSExpr, plot = '2d')
#> [[1]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> 
#> [[2]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> 
#> [[3]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

PcaPlot(dESeq2NormalizedDataFromTxLevelData$normalizedSExpr)
#> Warning in if (plot == "2d") {: the condition has length > 1 and only the first element will be used

#> Warning in if (plot == "2d") {: n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in if (plot == "2d") {: n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
#> Warning in normalizePath(f2): path[1]="webshot16db852124e0d.png": No such file
#> or directory
#> Warning in file(con, "rb"): cannot open file 'webshot16db852124e0d.png': No such
#> file or directory
#> Error in file(con, "rb"): cannot open the connection
PcaPlot(dESeq2NormalizedDataFromTxLevelData$normalizedSExpr, plot = '2d')
#> [[1]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> 
#> [[2]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> 
#> [[3]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

PcaPlot(tmmNormalizedData$normalizedSExpr)
#> Warning in if (plot == "2d") {: the condition has length > 1 and only the first element will be used

#> Warning in if (plot == "2d") {: n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in if (plot == "2d") {: n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
#> Warning in normalizePath(f2): path[1]="webshot16db85f2f97d4.png": No such file
#> or directory
#> Warning in file(con, "rb"): cannot open file 'webshot16db85f2f97d4.png': No such
#> file or directory
#> Error in file(con, "rb"): cannot open the connection
PcaPlot(tmmNormalizedData$normalizedSExpr, plot = '2d')
#> [[1]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> 
#> [[2]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> 
#> [[3]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

PcaPlot(tmmNormalizedDataFromTxLevelData$normalizedSExpr)
#> Warning in if (plot == "2d") {: the condition has length > 1 and only the first element will be used

#> Warning in if (plot == "2d") {: n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in if (plot == "2d") {: n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
#> Warning in normalizePath(f2): path[1]="webshot16db864f150c8.png": No such file
#> or directory
#> Warning in file(con, "rb"): cannot open file 'webshot16db864f150c8.png': No such
#> file or directory
#> Error in file(con, "rb"): cannot open the connection
PcaPlot(tmmNormalizedDataFromTxLevelData$normalizedSExpr, plot = '2d')
#> [[1]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> 
#> [[2]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
#> 
#> [[3]]
#> No trace type specified:
#>   Based on info supplied, a 'scatter' trace seems appropriate.
#>   Read more about this trace type -> https://plot.ly/r/reference/#scatter
#> No scatter mode specifed:
#>   Setting the mode to markers
#>   Read more about this attribute -> https://plot.ly/r/reference/#scatter-mode
#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors

#> Warning in RColorBrewer::brewer.pal(N, "Set2"): n too large, allowed maximum for palette Set2 is 8
#> Returning the palette you asked for with that many colors
```
