def excludelist(problem_type):
    # As of 20230426, the objective function of HS67 takes infinite time to be evaluated at some
    # points, e.g., [88.1351318; 12829.9219; 1.0e-5], maybe due to an infinite cycling.
    excludelist = ['HS76']
    # ARGLALE results in a segfault from slsqp when trying to project x0. We disable it for now
    excludelist += ['ARGLALE', 'ARGLBLE', 'ARGLCLE']

    if problem_type == 'unconstrained':
        # With the following excludelist, there is no unconstrained problem in CUTEst (as of 20230130) with
        # dimension between 51 and 100.
        excludelist += [
            'ARGTRIGLS',
            'BA-L1LS',
            'BA-L1SPLS',
            'BROWNAL',
            'CHNROSNB',
            'CHNRSNBM',
            'DIAMON2DLS',
            'DIAMON3DLS',
            'DMN15102LS',
            'DMN15103LS',
            'DMN15332LS',
            'DMN15333LS',
            'DMN37142LS',
            'DMN37143LS',
            'ERRINROS',
            'ERRINRSM',
            'HYDC20LS',
            'LUKSAN11LS',
            'LUKSAN12LS',
            'LUKSAN13LS',
            'LUKSAN14LS',
            'LUKSAN15LS',
            'LUKSAN16LS',
            'LUKSAN17LS',
            'LUKSAN22LS',
            'LUKSAN21LS',
            'LRCOVTYPE',
            'MANCINO',
            'QING',
            'SENSORS',
            'TOINTGOR',
            'TOINTPSP',
            'VARDIM',
        ]
    elif problem_type == 'bobyqa':
        # For the following problems, the classical bobyqa (single-precision) encounters SEGFAULT.
        excludelist += ['MGH17LS']
    elif problem_type == 'lincoa':
        excludelist += [
            'DALLASM', 
            'TARGUS',
        ]
    elif problem_type == 'cobyla':
        # The following problems were observed to take excessive time during tests GitHub Actions and
        # make the tests run overtime. Some of them may not be very time-consuming during a "plain"
        # test but become more challenging with some perturbations or variations. The excessive time may
        # be also a result of infinite cycling encountered by the classical version of cobyla.
        # The number following the problem is the time in seconds taken by cobyla in a "plain" test on
        # 20230130, which tested all linearly and nonlinearly constrained problems with at most 100
        # variables and 10000 constraints. Bound-constrained or unconstrained problems were not tested.
        excludelist += [
            'ACOPP30' ,
            'ACOPR30',
            'AIRPORT',      # 73
            'BATCH',        # 20
            'CHANDHEQ',     # 17
            'CHEBYQADNE',   # 546
            'CHNRSBNE',     # 18
            'CHNRSNBMNE',   # 32
            'CORE1',        # 64
            'CRESC100',
            'CRESC132',
            'CVXQP1',       # 54
            'DALLASS',      # 3 (it takes a long time on GitHub Actions)
            'DECONVBNE',
            'DECONVC',      # In a test on 20230328, the classical cobyla encountered infinite cycling.
            'DECONVNE',
            'DIAMON2D',     # 1415
            'DIAMON3D',     # 3703
            'DMN15102',     # 887
            'DMN15103',     # 3205
            'DMN15332',     # 838
            'DMN15333',     # 1441
            'DMN37142',     # 857
            'DMN37143',     # 2406
            'DUAL1',        # 73
            'DUAL2',        # 30
            'DUAL4',
            'ERRINRSMNE',   # 16
            'FBRAIN2',
            'FBRAIN2NE',
            'FBRAIN3',
            'FEEDLOC',
            'GOULDQP1',
            'HAIFAM',       # 173
            'HIMMELBI',     # 100
            'HIMMELBJ',
            'HYDCAR20',     # 648
            'HYDCAR6',
            'KISSING2',
            'LAKES',        # 65
            'LEVYMONE',     # 15
            'LHAIFAM',
            'LINSPANH',     # 3 (it takes a long time on GitHub Actions)
            'LUKSAN11',
            'LUKSAN12',     # 563
            'LUKSAN13',     # 508
            'LUKSAN14',     # 23
            'LUKSAN15',     # 19
            'LUKSAN16',     # 17
            'LUKSAN17',     # 25
            'LUKSAN21',     # 13
            'LUKSAN22',     # 19
            'MANCINONE',
            'MSS1',         # 39
            'OET5',
            'OET6',
            'OET7',
            'QINGNE',
            'QPCBLEND' ,
            'SPANHYD',      # 15
            'SWOPF',        # 10
            'TAX13322',     # 5
            'TAXR13322',    # 5
            'TRO4X4',       # 30
            'VANDERM1',     # 72
            'VANDERM2',     # 72
            'VANDERM3',     # 76
            'VESUVIO',
            'VESUVIOU',
        ]
    return excludelist