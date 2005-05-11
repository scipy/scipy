### This script is used to "manually" generate _na_cephesmodule.c

import sys
from numarray.codegenerator import UfuncModule, make_stub 
from numarray.numarrayext import NumarrayExtension

func_tab = [
(              'airy',                        'airy',       'dxdddd',     ['vxvvvv'],                    ['fxffff', 'dxdddd']),
(              'airy',                  'cairy_wrap',       'DxDDDD',     ['vxvvvv'],                    ['FxFFFF', 'DxDDDD']),
(             'airye',                'cairy_wrap_e',       'DxDDDD',     ['vxvvvv'],['fxffff', 'dxdddd', 'FxFFFF', 'DxDDDD']),
(              'bdtr',                        'bdtr',        'iidxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(             'bdtrc',                       'bdtrc',        'iidxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(             'bdtri',                       'bdtri',        'iidxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'bdtrik',                'cdfbin2_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'bdtrin',                'cdfbin3_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(               'bei',                    'bei_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'beip',                   'beip_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(               'ber',                    'ber_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'berp',                   'berp_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(        'besselpoly',                  'besselpoly',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(              'beta',                        'beta',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(           'betainc',                      'incbet',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(        'betaincinv',                       'incbi',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'betaln',                       'lbeta',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(             'btdtr',                       'btdtr',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'btdtri',                       'incbi',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(           'btdtria',                'cdfbet3_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(           'btdtrib',                'cdfbet4_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(              'cbrt',                        'cbrt',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'chdtr',                       'chdtr',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(            'chdtrc',                      'chdtrc',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(            'chdtri',                      'chdtri',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(           'chdtriv',                'cdfchi3_wrap',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(            'chndtr',                'cdfchn1_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(         'chndtridf',                'cdfchn3_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(         'chndtrinc',                'cdfchn4_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(          'chndtrix',                'cdfchn2_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(             'cosdg',                       'cosdg',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'cosm1',                       'cosm1',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'cotdg',                       'cotdg',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'dawsn',                       'dawsn',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(            'ellipe',                       'ellpe',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(         'ellipeinc',                       'ellie',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(            'ellipj',                       'ellpj',      'ddxdddd',    ['vvxvvvv'],                  ['ffxffff', 'ddxdddd']),
(            'ellipk',                       'ellpk',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(         'ellipkinc',                       'ellik',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(               'erf',                   'cerf_wrap',          'DxD',        ['vxf'],                          ['FxF', 'DxD']),
(               'erf',                         'erf',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'erfc',                        'erfc',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'exp1',                  'cexp1_wrap',          'DxD',        ['vxf'],                          ['FxF', 'DxD']),
(              'exp1',                   'exp1_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'exp10',                       'exp10',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'exp2',                        'exp2',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'expi',                   'expi_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'expm1',                       'expm1',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'expn',                        'expn',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
# duplicate (              'fdtr',                  'cdff1_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(              'fdtr',                        'fdtr',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(             'fdtrc',                       'fdtrc',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(             'fdtri',                       'fdtri',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(          'fdtridfd',                  'cdff4_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(          'fdtridfn',                  'cdff3_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
# undefined cdff2_wrap  (            'fdtrix',                  'cdff2_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(           'fresnel',                'cfresnl_wrap',         'DxDD',       ['vxvv'],                        ['FxFF', 'DxDD']),
(           'fresnel',                      'fresnl',         'dxdd',       ['vxvv'],                        ['fxff', 'dxdd']),
(             'gamma',                       'Gamma',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'gamma',                 'cgamma_wrap',          'DxD',        ['vxf'],                          ['FxF', 'DxD']),
(          'gammainc',                        'igam',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(         'gammaincc',                       'igamc',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(      'gammainccinv',                       'igami',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(           'gammaln',               'clngamma_wrap',          'DxD',        ['vxf'],                          ['FxF', 'DxD']),
(           'gammaln',                        'lgam',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'gdtr',                        'gdtr',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(             'gdtr2',                'cdfgam1_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(             'gdtrc',                       'gdtrc',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(             'gdtri',                       'gdtri',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'gdtria',                'cdfgam4_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'gdtrib',                'cdfgam3_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'gdtrix',                'cdfgam2_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(           'hankel1',                 'cbesh_wrap1',         'dDxD',       ['vvxf'],                        ['fFxF', 'dDxD']),
(          'hankel1e',               'cbesh_wrap1_e',         'dDxD',       ['vvxf'],                        ['fFxF', 'dDxD']),
(           'hankel2',                 'cbesh_wrap2',         'dDxD',       ['vvxf'],                        ['fFxF', 'dDxD']),
(          'hankel2e',               'cbesh_wrap2_e',         'dDxD',       ['vvxf'],                        ['fFxF', 'dDxD']),
(            'hyp1f1',                'chyp1f1_wrap',        'ddDxD',      ['vvvxf'],                      ['ffFxF', 'ddDxD']),
(            'hyp1f1',                      'hyperg',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'hyp1f2',                       'onef2',      'ddddxdd',    ['vvvvxfv'],                  ['ffffxff', 'ddddxdd']),
(            'hyp2f0',                      'hyp2f0',      'dddixdd',    ['vvvvxfv'],                  ['ffffxff', 'ddddxdd']),
(            'hyp2f1',                'chyp2f1_wrap',       'dddDxD',     ['vvvvxf'],                    ['fffFxF', 'dddDxD']),
(            'hyp2f1',                      'hyp2f1',       'ddddxd',     ['vvvvxf'],                    ['ffffxf', 'ddddxd']),
(            'hyp3f0',                     'threef0',      'ddddxdd',    ['vvvvxfv'],                  ['ffffxff', 'ddddxdd']),
(            'hyperu',                   'hypU_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(                'i0',                          'i0',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(               'i0e',                         'i0e',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(                'i1',                          'i1',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(               'i1e',                         'i1e',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(           'it2i0k0',                'it2i0k0_wrap',         'dxdd',       ['vxvv'],                        ['fxff', 'dxdd']),
(           'it2j0y0',                'it2j0y0_wrap',         'dxdd',       ['vxvv'],                        ['fxff', 'dxdd']),
(        'it2struve0',             'it2struve0_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(            'itairy',                 'itairy_wrap',       'dxdddd',     ['vxvvvv'],                    ['fxffff', 'dxdddd']),
(            'iti0k0',                'it1i0k0_wrap',         'dxdd',       ['vxvv'],                        ['fxff', 'dxdd']),
(            'itj0y0',                'it1j0y0_wrap',         'dxdd',       ['vxvv'],                        ['fxff', 'dxdd']),
(      'itmodstruve0',           'itmodstruve0_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(         'itstruve0',              'itstruve0_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(                'iv',                  'cbesi_wrap',         'dDxD',       ['vvxf'],                        ['fFxF', 'dDxD']),
(                'iv',                          'iv',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(               'ive',                'cbesi_wrap_e',         'dDxD',       ['vvxf'],        ['ffxf', 'ddxd', 'fFxF', 'dDxD']),
(                'j0',                          'j0',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(                'j1',                          'j1',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(                'jn',                          'jn',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(                'jv',                  'cbesj_wrap',         'dDxD',       ['vvxf'],                        ['fFxF', 'dDxD']),
(                'jv',                          'jv',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(               'jve',                'cbesj_wrap_e',         'dDxD',       ['vvxf'],        ['ffxf', 'ddxd', 'fFxF', 'dDxD']),
(                'k0',                          'k0',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(               'k0e',                         'k0e',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(                'k1',                          'k1',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(               'k1e',                         'k1e',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(               'kei',                    'kei_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'keip',                   'keip_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(            'kelvin',                 'kelvin_wrap',       'dxDDDD',     ['vxvvvv'],                    ['fxFFFF', 'dxDDDD']),
(               'ker',                    'ker_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'kerp',                   'kerp_wrap',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(                'kn',                          'kn',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(           'kolmogi',                     'kolmogi',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(        'kolmogorov',                  'kolmogorov',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(                'kv',                  'cbesk_wrap',         'dDxD',       ['vvxf'],        ['ffxf', 'ddxd', 'fFxF', 'dDxD']),
(               'kve',                'cbesk_wrap_e',         'dDxD',       ['vvxf'],        ['ffxf', 'ddxd', 'fFxF', 'dDxD']),
(             'log1p',                       'log1p',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(              'lpmv',                    'pmv_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(         'mathieu_a',                'cem_cva_wrap',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(         'mathieu_b',                'sem_cva_wrap',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(       'mathieu_cem',                    'cem_wrap',       'dddxdd',     ['vvvxvv'],                    ['fffxff', 'dddxdd']),
(   'mathieu_modcem1',                   'mcm1_wrap',       'dddxdd',     ['vvvxvv'],                    ['fffxff', 'dddxdd']),
(   'mathieu_modcem2',                   'mcm2_wrap',       'dddxdd',     ['vvvxvv'],                    ['fffxff', 'dddxdd']),
(   'mathieu_modsem1',                   'msm1_wrap',       'dddxdd',     ['vvvxvv'],                    ['fffxff', 'dddxdd']),
(   'mathieu_modsem2',                   'msm2_wrap',       'dddxdd',     ['vvvxvv'],                    ['fffxff', 'dddxdd']),
(       'mathieu_sem',                    'sem_wrap',       'dddxdd',     ['vvvxvv'],                    ['fffxff', 'dddxdd']),
(       'modfresnelm', 'modified_fresnel_minus_wrap',         'dxDD',       ['vxvv'],                        ['fxFF', 'dxDD']),
(       'modfresnelp',  'modified_fresnel_plus_wrap',         'dxDD',       ['vxvv'],                        ['fxFF', 'dxDD']),
(         'modstruve',              'modstruve_wrap',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(             'nbdtr',                       'nbdtr',        'iidxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'nbdtrc',                      'nbdtrc',        'iidxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'nbdtri',                      'nbdtri',        'iidxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(           'nbdtrik',                'cdfnbn2_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(           'nbdtrin',                'cdfnbn3_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'ncfdtr',                'cdffnc1_wrap',       'ddddxd',     ['vvvvxf'],                    ['ffffxf', 'ddddxd']),
(           'ncfdtri',                'cdffnc2_wrap',       'ddddxd',     ['vvvvxf'],                    ['ffffxf', 'ddddxd']),
(        'ncfdtridfd',                'cdffnc4_wrap',       'ddddxd',     ['vvvvxf'],                    ['ffffxf', 'ddddxd']),
(        'ncfdtridfn',                'cdffnc3_wrap',       'ddddxd',     ['vvvvxf'],                    ['ffffxf', 'ddddxd']),
(         'ncfdtrinc',                'cdffnc5_wrap',       'ddddxd',     ['vvvvxf'],                    ['ffffxf', 'ddddxd']),
(            'nctdtr',                'cdftnc1_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(         'nctdtridf',                'cdftnc3_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(         'nctdtrinc',                'cdftnc4_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(          'nctdtrit',                'cdftnc2_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(              'ndtr',                        'ndtr',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'ndtri',                       'ndtri',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(          'nrdtrimn',                'cdfnor3_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(          'nrdtrisd',                'cdfnor4_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(          'obl_ang1',      'oblate_aswfa_nocv_wrap',      'ddddxdd',    ['vvvvxvv'],                  ['ffffxff', 'ddddxdd']),
(       'obl_ang1_cv',           'oblate_aswfa_wrap',     'dddddxdd',   ['vvvvvxvv'],                ['fffffxff', 'dddddxdd']),
(            'obl_cv',            'oblate_segv_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(          'obl_rad1',    'oblate_radial1_nocv_wrap',      'ddddxdd',    ['vvvvxvv'],                  ['ffffxff', 'ddddxdd']),
(       'obl_rad1_cv',         'oblate_radial1_wrap',     'dddddxdd',   ['vvvvvxvv'],                ['fffffxff', 'dddddxdd']),
(          'obl_rad2',    'oblate_radial2_nocv_wrap',      'ddddxdd',    ['vvvvxvv'],                  ['ffffxff', 'ddddxdd']),
(       'obl_rad2_cv',         'oblate_radial2_wrap',     'dddddxdd',   ['vvvvvxvv'],                ['fffffxff', 'dddddxdd']),
(              'pbdv',                   'pbdv_wrap',        'ddxdd',      ['vvxvv'],                      ['ffxff', 'ddxdd']),
(              'pbvv',                   'pbvv_wrap',        'ddxdd',      ['vvxvv'],                      ['ffxff', 'ddxdd']),
(              'pbwa',                   'pbwa_wrap',        'ddxdd',      ['vvxvv'],                      ['ffxff', 'ddxdd']),
(              'pdtr',                        'pdtr',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(             'pdtrc',                       'pdtrc',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(             'pdtri',                       'pdtri',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(            'pdtrik',                'cdfpoi2_wrap',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(          'pro_ang1',     'prolate_aswfa_nocv_wrap',      'ddddxdd',    ['vvvvxvv'],                  ['ffffxff', 'ddddxdd']),
(       'pro_ang1_cv',          'prolate_aswfa_wrap',     'dddddxdd',   ['vvvvvxvv'],                ['fffffxff', 'dddddxdd']),
(            'pro_cv',           'prolate_segv_wrap',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(          'pro_rad1',   'prolate_radial1_nocv_wrap',      'ddddxdd',    ['vvvvxvv'],                  ['ffffxff', 'ddddxdd']),
(       'pro_rad1_cv',        'prolate_radial1_wrap',     'dddddxdd',   ['vvvvvxvv'],                ['fffffxff', 'dddddxdd']),
(          'pro_rad2',   'prolate_radial2_nocv_wrap',      'ddddxdd',    ['vvvvxvv'],                  ['ffffxff', 'ddddxdd']),
(       'pro_rad2_cv',        'prolate_radial2_wrap',     'dddddxdd',   ['vvvvvxvv'],                ['fffffxff', 'dddddxdd']),
(               'psi',                   'cpsi_wrap',          'DxD',        ['vxf'],                          ['FxF', 'DxD']),
(               'psi',                         'psi',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(            'radian',                      'radian',        'dddxd',      ['vvvxf'],                      ['fffxf', 'dddxd']),
(            'rgamma',                'crgamma_wrap',          'DxD',        ['vxf'],                          ['FxF', 'DxD']),
(            'rgamma',                      'rgamma',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'round',                       'round',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(            'shichi',                      'shichi',         'dxdd',       ['vxvv'],                        ['fxff', 'dxdd']),
(              'sici',                        'sici',         'dxdd',       ['vxvv'],                        ['fxff', 'dxdd']),
(             'sindg',                       'sindg',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(           'smirnov',                     'smirnov',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(          'smirnovi',                    'smirnovi',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(            'spence',                      'spence',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(             'stdtr',                  'cdft1_wrap',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
# (             'stdtr',                       'stdtr',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(            'stdtri',                      'stdtri',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(          'stdtridf',                  'cdft3_wrap',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(           'stdtrit',                  'cdft2_wrap',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(            'struve',                      'struve',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(             'tandg',                       'tandg',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(           'tklmbda',              'tukeylambdacdf',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(              'wofz',                  'cwofz_wrap',          'DxD',        ['vxf'],                          ['FxF', 'DxD']),
(                'y0',                          'y0',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(                'y1',                          'y1',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
(                'yn',                          'yn',         'idxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(                'yv',                  'cbesy_wrap',         'dDxD',       ['vvxf'],                        ['fFxF', 'dDxD']),
(                'yv',                          'yv',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(               'yve',                'cbesy_wrap_e',         'dDxD',       ['vvxf'],        ['ffxf', 'ddxd', 'fFxF', 'dDxD']),
(              'zeta',                        'zeta',         'ddxd',       ['vvxf'],                        ['ffxf', 'ddxd']),
(             'zetac',                       'zetac',          'dxd',        ['vxf'],                          ['fxf', 'dxd']),
]


python_code = '''

_scipy_special_errprint = 0

def errprint(val=None):
    global _scipy_special_errprint
    old_val = _scipy_special_errprint
    if val is not None:
        _scipy_special_errprint = (val != 0)
    return old_val
'''

# Added "fake" headers due to naming conflicts between Py_complex
# (really used by cephes wrappers) and Complex{32,64} (used by numarray).

fake_header_code = '''#include "cephes.h"

double besselpoly(double a, double lambda, double nu);

Complex64 cbesi_wrap( double v, Complex64 z);
Complex64 cbesi_wrap_e( double v, Complex64 z);
Complex64 cbesj_wrap( double v, Complex64 z);
Complex64 cbesj_wrap_e( double v, Complex64 z);
Complex64 cbesy_wrap( double v, Complex64 z);
Complex64 cbesy_wrap_e( double v, Complex64 z);
Complex64 cbesk_wrap( double v, Complex64 z);
Complex64 cbesk_wrap_e( double v, Complex64 z);  
Complex64 cbesh_wrap1( double v, Complex64 z);
Complex64 cbesh_wrap1_e( double v, Complex64 z);  
Complex64 cbesh_wrap2( double v, Complex64 z);
Complex64 cbesh_wrap2_e( double v, Complex64 z);

extern double cdfbet3_wrap(double p, double x, double b);
extern double cdfbet4_wrap(double p, double x, double a);

extern double cdfbin2_wrap(double p, double xn, double pr);
extern double cdfbin3_wrap(double p, double s, double pr);

extern double cdfchi3_wrap(double p, double x);

extern double cdfchn1_wrap(double x, double df, double nc);
extern double cdfchn2_wrap(double p, double df, double nc);
extern double cdfchn3_wrap(double p, double x, double nc);
extern double cdfchn4_wrap(double p, double x, double df);

extern double cdff3_wrap(double p, double f, double dfd);
extern double cdff4_wrap(double p, double f, double dfn);

extern double cdffnc1_wrap(double f, double dfn, double dfd, double nc);
extern double cdffnc2_wrap(double p, double dfn, double dfd, double nc);
extern double cdffnc3_wrap(double p, double f, double dfd, double nc);
extern double cdffnc4_wrap(double p, double f, double dfn, double nc);
extern double cdffnc5_wrap(double p, double f, double dfn, double dfd);

extern double cdfgam1_wrap(double p, double x, double scl);
extern double cdfgam2_wrap(double p, double x, double shp);
extern double cdfgam3_wrap(double p, double x, double scl);
extern double cdfgam4_wrap(double p, double x, double shp);

extern double cdfnbn2_wrap(double p, double xn, double pr);
extern double cdfnbn3_wrap(double p, double s, double pr);

extern double cdfnor3_wrap(double p, double x, double std);
extern double cdfnor4_wrap(double p, double x, double mn);

extern double cdfpoi2_wrap(double p, double xlam);

extern double cdft1_wrap(double p, double t);
extern double cdft2_wrap(double p, double t);
extern double cdft3_wrap(double p, double t);

extern double cdftnc1_wrap(double df, double nc, double t);
extern double cdftnc2_wrap(double df, double nc, double p);
extern double cdftnc3_wrap(double p, double nc, double t);
extern double cdftnc4_wrap(double df, double p, double t);

extern double tukeylambdacdf(double x, double lambda);

Complex64 cgamma_wrap( Complex64 z);
Complex64 clngamma_wrap( Complex64 z);
Complex64 cpsi_wrap( Complex64 z);
Complex64 crgamma_wrap( Complex64 z);
Complex64 chyp2f1_wrap( double a, double b, double c, Complex64 z);
Complex64 chyp1f1_wrap( double a, double b, Complex64 z);
double hypU_wrap(double a, double b, double x);
double exp1_wrap(double x);
double expi_wrap(double x);
Complex64 cexp1_wrap( Complex64 z);
Complex64 cerf_wrap( Complex64 z);
int itairy_wrap(double x, double *apt, double *bpt, double *ant, double *bnt);

double struve_wrap(double v, double x);
double itstruve0_wrap(double x);
double it2struve0_wrap(double x);

double modstruve_wrap(double v, double x);
double itmodstruve0_wrap(double x);

double ber_wrap(double x);
double bei_wrap(double x);
double ker_wrap(double x);
double kei_wrap(double x);
double berp_wrap(double x);
double beip_wrap(double x);
double kerp_wrap(double x);
double keip_wrap(double x);

int kelvin_wrap(double x, Complex64 *Be, Complex64 *Ke, Complex64 *Bep, Complex64 *Kep);

int it1j0y0_wrap(double x, double *, double *);
int it2j0y0_wrap(double x, double *, double *);
int it1i0k0_wrap(double x, double *, double *);
int it2i0k0_wrap(double x, double *, double *);

int cfresnl_wrap(Complex64 x, Complex64 *sf, Complex64 *cf);
double cem_cva_wrap(double m, double q);
double sem_cva_wrap(double m, double q);
int cem_wrap(double m, double q, double x, double *csf, double *csd);
int sem_wrap(double m, double q, double x, double *csf, double *csd);
int mcm1_wrap(double m, double q, double x, double *f1r, double *d1r);
int msm1_wrap(double m, double q, double x, double *f1r, double *d1r);
int mcm2_wrap(double m, double q, double x, double *f2r, double *d2r);
int msm2_wrap(double m, double q, double x, double *f2r, double *d2r);
double pmv_wrap(double, double, double);
int pbwa_wrap(double, double, double *, double *);
int pbdv_wrap(double, double, double *, double *);
int pbvv_wrap(double, double, double *, double *);

int prolate_aswfa_wrap(double, double, double, double, double, double *, double *);
int prolate_radial1_wrap(double, double, double, double, double, double *, double *);
int prolate_radial2_wrap(double, double, double, double, double, double *, double *);

/*
int oblate_aswfa_wrap(double, double, double, double, double, double *, double *);
int oblate_radial1_wrap(double, double, double, double, double, double *, double *);
int oblate_radial2_wrap(double, double, double, double, double, double *, double *);
double prolate_aswfa_nocv_wrap(double, double, double, double, double *);
double prolate_radial1_nocv_wrap(double, double, double, double, double *);
double prolate_radial2_nocv_wrap(double, double, double, double, double *);
double oblate_aswfa_nocv_wrap(double, double, double, double, double *);
double oblate_radial1_nocv_wrap(double, double, double, double, double *);
double oblate_radial2_nocv_wrap(double, double, double, double, double *);
*/

double prolate_segv_wrap(double, double, double);
double oblate_segv_wrap(double, double, double);


int modified_fresnel_plus_wrap(double x, Complex64 *F, Complex64 *K);
int modified_fresnel_minus_wrap(double x, Complex64 *F, Complex64 *K);

extern Complex64 cwofz_wrap(Complex64 z);

'''


def gencode():
    m = UfuncModule('_cephes')

    m.add_code( fake_header_code ) 
    
    for r in func_tab:
        ufunc, c_func, c_sig, forms, ufsigs = r
        print "adding",ufunc
        m.add_nary_ufunc( ufunc_name=ufunc,
                          c_function=c_func,
                          forms=forms,
                          signatures=ufsigs,
                          c_signature=c_sig)
        
    m.generate('_na_cephesmodule.c')

    make_stub('na_cephes', '_cephes',
              add_code = python_code)  

if __name__ == "__main__":
    gencode()

