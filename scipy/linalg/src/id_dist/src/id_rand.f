c       this file contains the following user-callable routines:
c
c
c       routine id_frand generates pseudorandom numbers
c       drawn uniformly from [0,1]. id_frand is more
c       efficient that id_srand, but cannot generate
c       fewer than 55 pseudorandom numbers per call.
c
c       routine id_srand generates pseudorandom numbers
c       drawn uniformly from [0,1]. id_srand is less
c       efficient that id_frand, but can generate
c       fewer than 55 pseudorandom numbers per call.
c
c       entry id_frandi initializes the seed values
c       for routine id_frand.
c
c       entry id_srandi initializes the seed values
c       for routine id_srand.
c
c       entry id_frando initializes the seed values
c       for routine id_frand to their original values.
c
c       entry id_srando initializes the seed values
c       for routine id_srand to their original values.
c
c       routine id_randperm generates a uniformly random permutation.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine id_frand(n,r)
c
c       generates n pseudorandom numbers drawn uniformly from [0,1],
c       via a very efficient lagged Fibonnaci method.
c       Unlike routine id_srand, the present routine requires that
c       n be at least 55.
c
c       input:
c       n -- number of pseudorandom numbers to generate
c
c       output:
c       r -- array of pseudorandom numbers
c
c       _N.B._: n must be at least 55.
c
c       reference:
c       Press, Teukolsky, Vetterling, Flannery, "Numerical Recipes,"
c            3rd edition, Cambridge University Press, 2007,
c            Section 7.1.5.
c
        implicit none
        integer n,k
        real*8 r(n),s(55),t(55),s0(55),x
        save
c
        data s/
     1  0.2793574644042651d0, 0.1882566493961346d0,
     2  0.5202478134503912d0, 0.7568505373052146d0,
     3  0.5682465992936152d0, 0.5153148754383294d0,
     4  0.7806554095454596d0, 1.982474428974643d-2,
     5  0.2520464262278498d0, 0.6423784715775962d0,
     6  0.5802024387972178d0, 0.3784471040388249d0,
     7  7.839919528229308d-2, 0.6334519212594525d0,
     8  3.387627157788001d-2, 0.1709066283884670d0,
     9  0.4801610983518325d0, 0.8983424668099422d0,
     *  5.358948687598758d-2, 0.1265377231771848d0,
     1  0.8979988627693677d0, 0.6470084038238917d0,
     2  0.3031709395541237d0, 0.6674702804438126d0,
     3  0.6318240977112699d0, 0.2235229633873050d0,
     4  0.2784629939177633d0, 0.2365462014457445d0,
     5  0.7226213454977284d0, 0.8986523045307989d0,
     6  0.5488233229247885d0, 0.3924605412141200d0,
     7  0.6288356378374988d0, 0.6370664115760445d0,
     8  0.5925600062791174d0, 0.4322113919396362d0,
     9  0.9766098520360393d0, 0.5168619893947437d0,
     *  0.6799970440779681d0, 0.4196004604766881d0,
     1  0.2324473089903044d0, 0.1439046416143282d0,
     2  0.4670307948601256d0, 0.7076498261128343d0,
     3  0.9458030397562582d0, 0.4557892460080424d0,
     4  0.3905930854589403d0, 0.3361770064397268d0,
     5  0.8303274937900278d0, 0.3041110304032945d0,
     6  0.5752684022049654d0, 7.985703137991175d-2,
     7  0.5522643936454465d0, 1.956754937251801d-2,
     8  0.9920272858340107d0/
c
        data s0/
     1  0.2793574644042651d0, 0.1882566493961346d0,
     2  0.5202478134503912d0, 0.7568505373052146d0,
     3  0.5682465992936152d0, 0.5153148754383294d0,
     4  0.7806554095454596d0, 1.982474428974643d-2,
     5  0.2520464262278498d0, 0.6423784715775962d0,
     6  0.5802024387972178d0, 0.3784471040388249d0,
     7  7.839919528229308d-2, 0.6334519212594525d0,
     8  3.387627157788001d-2, 0.1709066283884670d0,
     9  0.4801610983518325d0, 0.8983424668099422d0,
     *  5.358948687598758d-2, 0.1265377231771848d0,
     1  0.8979988627693677d0, 0.6470084038238917d0,
     2  0.3031709395541237d0, 0.6674702804438126d0,
     3  0.6318240977112699d0, 0.2235229633873050d0,
     4  0.2784629939177633d0, 0.2365462014457445d0,
     5  0.7226213454977284d0, 0.8986523045307989d0,
     6  0.5488233229247885d0, 0.3924605412141200d0,
     7  0.6288356378374988d0, 0.6370664115760445d0,
     8  0.5925600062791174d0, 0.4322113919396362d0,
     9  0.9766098520360393d0, 0.5168619893947437d0,
     *  0.6799970440779681d0, 0.4196004604766881d0,
     1  0.2324473089903044d0, 0.1439046416143282d0,
     2  0.4670307948601256d0, 0.7076498261128343d0,
     3  0.9458030397562582d0, 0.4557892460080424d0,
     4  0.3905930854589403d0, 0.3361770064397268d0,
     5  0.8303274937900278d0, 0.3041110304032945d0,
     6  0.5752684022049654d0, 7.985703137991175d-2,
     7  0.5522643936454465d0, 1.956754937251801d-2,
     8  0.9920272858340107d0/
c
c
        do k = 1,24
c
          x = s(k+31)-s(k)
          if(x .lt. 0) x = x+1
          r(k) = x
c
        enddo ! k
c
c
        do k = 25,55
c
          x = r(k-24)-s(k)
          if(x .lt. 0) x = x+1
          r(k) = x
c
        enddo ! k
c
c
        do k = 56,n
c
          x = r(k-24)-r(k-55)
          if(x .lt. 0) x = x+1
          r(k) = x
c
        enddo ! k
c
c
        do k = 1,55
          s(k) = r(n-55+k)
        enddo ! k
c
c
        return
c
c
c
        entry id_frandi(t)
c
c       initializes the seed values in s
c       (any appropriately random numbers will do).
c
c       input:
c       t -- values to copy into s
c
        do k = 1,55
          s(k) = t(k)
        enddo ! k
c
        return
c
c
c
        entry id_frando()
c
c       initializes the seed values in s to their original values.
c
        do k = 1,55
          s(k) = s0(k)
        enddo ! k
c
        return
        end
c
c
c
c
        subroutine id_srand(n,r)
c
c       generates n pseudorandom numbers drawn uniformly from [0,1],
c       via a very efficient lagged Fibonnaci method.
c       Unlike routine id_frand, the present routine does not requires
c       that n be at least 55.
c
c       input:
c       n -- number of pseudorandom numbers to generate
c
c       output:
c       r -- array of pseudorandom numbers
c
c       reference:
c       Press, Teukolsky, Vetterling, Flannery, "Numerical Recipes,"
c            3rd edition, Cambridge University Press, 2007,
c            Section 7.1.5.
c
        implicit none
        integer n,k,l,m
        real*8 s(55),r(n),s0(55),t(55),x
        save
c
        data l/55/,m/24/
c
        data s/
     1  0.8966049453474352d0, 0.7789471911260157d0,
     2  0.6071529762908476d0, 0.8287077988663865d0,
     3  0.8249336255502409d0, 0.5735259423199479d0,
     4  0.2436346323812991d0, 0.2656149927259701d0,
     5  0.6594784809929011d0, 0.3432392503145575d0,
     6  0.5051287353012308d0, 0.1444493249757482d0,
     7  0.7643753221285416d0, 0.4843422506977382d0,
     8  0.4427513254774826d0, 0.2965991475108561d0,
     9  0.2650513544474467d0, 2.768759325778929d-2,
     *  0.6106305243078063d0, 0.4246918885003141d0,
     1  0.2863757386932874d0, 0.6211983878375777d0,
     2  0.7534336463880467d0, 0.7471458603576737d0,
     3  0.2017455446928328d0, 0.9334235874832779d0,
     4  0.6343440435422822d0, 0.8819824804812527d0,
     5  1.994761401222460d-2, 0.7023693520374801d0,
     6  0.6010088924817263d0, 6.498095955562046d-2,
     7  0.3090915456102685d0, 0.3014924769096677d0,
     8  0.5820726822705102d0, 0.3630527222866207d0,
     9  0.3787166916242271d0, 0.3932772088505305d0,
     *  0.5570720335382000d0, 0.9712062146993835d0,
     1  0.1338293907964648d0, 0.1857441593107195d0,
     2  0.9102503893692572d0, 0.2623337538798778d0,
     3  0.3542828591321135d0, 2.246286032456513d-2,
     4  0.7935703170405717d0, 6.051464729640567d-2,
     5  0.7271929955172147d0, 1.968513010678739d-3,
     6  0.4914223624495486d0, 0.8730023176789450d0,
     7  0.9639777091743168d0, 0.1084256187532446d0,
     8  0.8539399636754000d0/
c
        data s0/
     1  0.8966049453474352d0, 0.7789471911260157d0,
     2  0.6071529762908476d0, 0.8287077988663865d0,
     3  0.8249336255502409d0, 0.5735259423199479d0,
     4  0.2436346323812991d0, 0.2656149927259701d0,
     5  0.6594784809929011d0, 0.3432392503145575d0,
     6  0.5051287353012308d0, 0.1444493249757482d0,
     7  0.7643753221285416d0, 0.4843422506977382d0,
     8  0.4427513254774826d0, 0.2965991475108561d0,
     9  0.2650513544474467d0, 2.768759325778929d-2,
     *  0.6106305243078063d0, 0.4246918885003141d0,
     1  0.2863757386932874d0, 0.6211983878375777d0,
     2  0.7534336463880467d0, 0.7471458603576737d0,
     3  0.2017455446928328d0, 0.9334235874832779d0,
     4  0.6343440435422822d0, 0.8819824804812527d0,
     5  1.994761401222460d-2, 0.7023693520374801d0,
     6  0.6010088924817263d0, 6.498095955562046d-2,
     7  0.3090915456102685d0, 0.3014924769096677d0,
     8  0.5820726822705102d0, 0.3630527222866207d0,
     9  0.3787166916242271d0, 0.3932772088505305d0,
     *  0.5570720335382000d0, 0.9712062146993835d0,
     1  0.1338293907964648d0, 0.1857441593107195d0,
     2  0.9102503893692572d0, 0.2623337538798778d0,
     3  0.3542828591321135d0, 2.246286032456513d-2,
     4  0.7935703170405717d0, 6.051464729640567d-2,
     5  0.7271929955172147d0, 1.968513010678739d-3,
     6  0.4914223624495486d0, 0.8730023176789450d0,
     7  0.9639777091743168d0, 0.1084256187532446d0,
     8  0.8539399636754000d0/
c
c
        do k = 1,n
c
c         Run one step of the recurrence.
c
          x = s(m)-s(l)
          if(x .lt. 0) x = x+1
          s(l) = x
          r(k) = x
c
c         Decrement l and m.
c
          l = l-1
          m = m-1
c
c         Circle back to the end if required.
c
          if(l .eq. 0) l = 55
          if(m .eq. 0) m = 55
c
        enddo ! k
c
c
        return
c
c
c
        entry id_srandi(t)
c
c       initializes the seed values in s
c       (any appropriately random numbers will do).
c
c       input:
c       t -- values to copy into s
c
        do k = 1,55
          s(k) = t(k)
        enddo ! k
c
        l = 55
        m = 24
c
        return
c
c
c
        entry id_srando()
c
c       initializes the seed values in s to their original values.
c
        do k = 1,55
          s(k) = s0(k)
        enddo ! k
c
        l = 55
        m = 24
c
        return
        end
c
c
c
c
        subroutine id_randperm(n,ind)
c
c       draws a permutation ind uniformly at random from the group
c       of all permutations of n objects.
c
c       input:
c       n -- length of ind
c
c       output:
c       ind -- random permutation of length n
c
        implicit none
        integer n,ind(n),m,j,iswap
        real*8 r
c
c
c       Initialize ind.
c
        do j = 1,n
          ind(j) = j
        enddo ! j
c
c
c       Shuffle ind via the Fisher-Yates (Knuth/Durstenfeld) algorithm.
c
        do m = n,2,-1
c
c         Draw an integer uniformly at random from 1, 2, ..., m.
c
          call id_srand(1,r)
          j = m*r+1
c
c         Uncomment the following line if r could equal 1:
c         if(j .eq. m+1) j = m
c
c         Swap ind(j) and ind(m).
c
          iswap = ind(j)
          ind(j) = ind(m)
          ind(m) = iswap
c
        enddo ! m
c
c
        return
        end
