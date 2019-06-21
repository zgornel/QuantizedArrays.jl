###################################################################################
#    'dko'                                         ,,                        :NM. #
#  'Ol..'lO'                                 d,    od                         .M. #
# ,W'     ,W,              ..         ..     Xl                   ...      .. .M. #
# K0       K0 ,OX  'kW   ox,,dk   cNklckX.  cNOcc ;OX  lxcclWK  '0c,l0.  .Kx:cdM. #
# Wk       ON  cX   :W   .    Xc   Kk   Xo   Xc    cX  ..  OO  .W:   lN  Xx   ,M. #
# X0       KK  cX   :W   cOollNl   Kl   0o   Xc    cX    .Kx   ;MdcccoO..Mc   .M. #
# :W.     .M;  cN   cW  'M,   Kl   Kl   0o   Xc    cX   .Nl    'M,    '  Wo   'M. #
#  ;0:. .:0;   .Wl.;OW.  Xx. ;Wd. .Xd. .Kx.  0d d,.dN. ,Wd...x' l0. .do  :N;..xM:.#
#    ':oNc      .:c..,,   ;c:.',' ,,,. ,,,'  .:c, ',,, .,,,,,,   .:c:.    .:c, ,,.#
#       .kK,                                                                      #
#                                                                                 #
#                        ....      .....               ...   ...                  #
#          ,::::::;     ,MMMM,     XMMMx     :::::::. ;MMo   kMM.    '::;::,      #
#          XMMMMMMX     ,MMMM,     XMMMx    .MMMMMMMd ;MMo   kMM.    0MMWMM0      #
#          XMMMMMMX  kOOo''''  cOOO,'''.    .MMMMMMMd  '':OOO;''     0MMWMM0      #
#          XMMMMMMX  WMMx      xMMM         .MMMMMMMd    :MMM'       0MMX         #
#          XMMMMMMX  WMMx      xMMM         .MMMMMMMd    :MMM'       0MMX         #
#       OWWWMMMMMMX  WMMx      xMMM      'WWNMMMMMMMd    :MMM'    xWWWMMX         #
#       0MMWMMMMMMX  WMMx      xMMM      'MMWMMMMMMMd    :MMM'    kMMWMMX         #
#       0MMWMMMMMMX  WMMx      xMMM      'MMWMMMMMMMd    :MMM'    kMMWMMX         #
#                                                                                 #
###################################################################################

# QuantizedArrays.jl - Array quantization and compression, written at 0x0α Research
#                      by Corneliu Cofaru, 2019

__precompile__(true)

"""
QuantizedArrays.jl - array quantization and compression.
"""
module QuantizedArrays

using LinearAlgebra,
      InteractiveUtils,
      StatsBase,
      Distances,
      Clustering

export QuantizedArray,
       QuantizedVector,
       QuantizedMatrix,
       ArrayQuantizer,
       CodeBook,
       quantize,
       build_quantizer,
       build_codebooks

include("defaults.jl")
include("utils.jl")
include("codebook.jl")
include("codebook_algorithms.jl")
include("quantizer.jl")
include("interface.jl")

end # module
