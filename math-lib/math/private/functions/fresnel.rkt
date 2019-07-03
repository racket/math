#lang typed/racket/base

#|
Fresnel Integrals:
 S(x)=∫sin(π/2t²) |0->x
 C(x)=∫cos(π/2t²) |0->x

floating point implementation:
  adaptation of the algorithm in the fresnl.c file from Cephes,
  but with the limit for the rational powerseries lowered to x<1.5625 (coming from 2.5625)
  above x>1.6 the error becomes too big (~1e-8 x=2.5)
Bigfloat implementation:
  adaptation of powerseries as shown on wikipedia

would like to extend to the complex plane
S(xi)=-S(x)i
C(xi)=C(x)i

S((1+i)/√2 × x)=(-1+i)/√2 × S(x)
S((-1+i)/√2 × x)=(1+i)/√2 × S(x)

checking also possible via erf((1+i)/√2 × z)

TODO:
put tests in math-test/math/tests/function-tests.rkt
update documentation
|#

(require "../../base.rkt"
         "../../flonum.rkt")

(provide flFresnel-S Fresnel-S flFresnel-C Fresnel-C)

;------------------------------
;polinomials for the rational assymptotic aproximation for S & C
(define fn (make-flpolyfun ( 3.76329711269987889006E-20
                             1.34283276233062758925E-16
                             1.72010743268161828879E-13
                             1.02304514164907233465E-10
                             3.05568983790257605827E-8
                             4.63613749287867322088E-6
                             3.45017939782574027900E-4
                             1.15220955073585758835E-2
                             1.43407919780758885261E-1
                             4.21543555043677546506E-1)))
(define fd (make-flpolyfun ( 1.25443237090011264384E-20
                             4.52001434074129701496E-17
                             5.88754533621578410010E-14
                             3.60140029589371370404E-11
                             1.12699224763999035261E-8
                             1.84627567348930545870E-6
                             1.55934409164153020873E-4
                             6.44051526508858611005E-3
                             1.16888925859191382142E-1
                             7.51586398353378947175E-1
                             1.0)))
(define gn (make-flpolyfun ( 1.86958710162783235106E-22
                             8.36354435630677421531E-19
                             1.37555460633261799868E-15
                             1.08268041139020870318E-12
                             4.45344415861750144738E-10
                             9.82852443688422223854E-8
                             1.15138826111884280931E-5
                             6.84079380915393090172E-4
                             1.87648584092575249293E-2
                             1.97102833525523411709E-1
                             5.04442073643383265887E-1)))
(define gd (make-flpolyfun ( 1.86958710162783236342E-22
                             8.39158816283118707363E-19
                             1.38796531259578871258E-15
                             1.10273215066240270757E-12
                             4.60680728146520428211E-10
                             1.04314589657571990585E-7
                             1.27545075667729118702E-5
                             8.14679107184306179049E-4
                             2.53603741420338795122E-2
                             3.37748989120019970451E-1
                             1.47495759925128324529E0
                             1.0)))

;------------------------------
;polinomials for the rational powerseries aproximation for S
(define sn (make-flpolyfun ( 3.18016297876567817986E11
                            -4.42979518059697779103E10
                             2.54890880573376359104E9
                            -6.29741486205862506537E7
                             7.08840045257738576863E5
                            -2.99181919401019853726E3)))
(define sd (make-flpolyfun ( 6.07366389490084639049E11
                             2.24411795645340920940E10
                             4.19320245898111231129E8
                             5.17343888770096400730E6
                             4.55847810806532581675E4
                             2.81376268889994315696E2
                             1.0)))
(define (flFresnel-S [x : Float]) : Float
  (cond
    [(fl< x 0.0) (fl* -1.0 (flFresnel-S (fl* -1.0 x)))]
    [(fl= x 0.0) 0.0]
    [(fl< 36974.0 x) .5]
    [(fl< x 1.5625)
     (define X2 (fl* x x))
     (define X4 (fl* X2 X2))
     (fl* (fl* X2 x) (fl/ (sn X4)(sd X4)))]
    [else
     (define t (fl* pi (fl* x x)))
     (define t/2 (fl/ t 2.0))
     (define U (fl/ 1.0 (fl* t t)))
     (define f (fl- 1.0 (fl* U (fl/ (fn U)(fd U)))))
     (define g (fl/ (fl/ (gn U)(gd U)) t))
     (fl- 0.5 (fl/ (fl+ (fl* f (flcos t/2))(fl* g (flsin t/2)))
                   (fl* pi x)))]))

;------------------------------
;polinomials for the rational powerseries aproximation for C
(define cn (make-flpolyfun ( 9.99999999999999998822E-1
                            -2.05525900955013891793E-1
                             1.88843319396703850064E-2
                            -6.45191435683965050962E-4
                             9.50428062829859605134E-6
                            -4.98843114573573548651E-8)))
(define cd (make-flpolyfun ( 1.00000000000000000118E0
                             4.12142090722199792936E-2
                             8.68029542941784300606E-4
                             1.22262789024179030997E-5
                             1.25001862479598821474E-7
                             9.15439215774657478799E-10
                             3.99982968972495980367E-12)))
(define (flFresnel-C [x : Float]) : Float
  (cond
    [(fl< x 0.0) (fl* -1.0 (flFresnel-C (fl* -1.0 x)))]
    [(fl= x 0.0) 0.0]
    [(fl< 36974.0 x) .5]
    [(fl< x 1.5625)
     (define X2 (fl* x x))
     (define X4 (fl* X2 X2))
     (fl* x (fl/ (cn X4)(cd X4)))]
    [else
     (define t (fl* pi (fl* x x)))
     (define t/2 (fl/ t 2.0))
     (define U (fl/ 1.0 (fl* t t)))
     (define f (fl- 1.0 (fl* U (fl/ (fn U)(fd U)))))
     (define g (fl/ (fl/ (gn U)(gd U)) t))
     (fl+ 0.5 (fl/ (fl- (fl* f (flsin t/2))(fl* g (flcos t/2)))
                   (fl* pi x)))]))

;------------------------------
(define (Fresnel-S [x : Real]) : Real (flFresnel-S (fl x)))
(define (Fresnel-C [x : Real]) : Real (flFresnel-C (fl x)))

;------------------------------
(module* bfFresnel #f
  (require math/bigfloat)
  (provide bfFresnel-S bfFresnel-C)
  ;this implementation is potentially exact
  ;but for large x (> 5!) it is really slow and needs a lot of bits in bf-precision!
  (define (bfFresnel-S [x : Bigfloat])
    (define X (bf* x (bfsqrt (bf/ pi.bf (bf 2)))))
    (define X4 (bfexpt X (bf 4)))
    (define prsn (bfexpt (bf 1/2) (bf (bf-precision))))
    (define-values (s l)
       (for/fold : (Values Bigfloat Bigfloat)
         ([s : Bigfloat (bf/ (bfexpt X (bf 3)) (bf 3))]
          [l : Bigfloat (bf/ (bfexpt X (bf 3)) (bf 3))])
         ([n (in-naturals 1)]
          #:break (bf< (bfabs (bf/ l s)) prsn))
         (define l+ (bf* (bf -1) X4 (bf* l (bf (/ (- (* 4 n) 1)
                                                  (* 2 n)(+ (* 2 n) 1)(+ (* 4 n) 3))))))
         (values (bf+ s l+) l+)))
     (bf* (bfsqrt (bf/ (bf 2) pi.bf)) s))
  (define (bfFresnel-C [x : Bigfloat])
    (define X (bf* x (bfsqrt (bf/ pi.bf (bf 2)))))
    (define X4 (bfexpt X (bf 4)))
    (define prsn (bfexpt (bf 1/2) (bf (bf-precision))))
    (define-values (s l)
       (for/fold : (Values Bigfloat Bigfloat)
         ([s : Bigfloat (bf/ (bfexpt X (bf 1)) (bf 1))]
          [l : Bigfloat (bf/ (bfexpt X (bf 1)) (bf 1))])
         ([n (in-naturals 1)]
          #:break (bf< (bfabs (bf/ l s)) prsn))
         (define l+ (bf* (bf -1) X4 (bf* l (bf (/ (- (* 4 n) 3)
                                                  (* 2 n)(- (* 2 n) 1)(+ (* 4 n) 1))))))
         (values (bf+ s l+) l+)))
     (bf* (bfsqrt (bf/ (bf 2) pi.bf)) s))
  )

#;(module* test racket/base
  (require math/bigfloat
           rackunit)
  (require (submod "..")
           (submod ".." bfFresnel))
  
  (define (test a)
    (check-= (Fresnel-S a) (bigfloat->flonum (bfFresnel-S (bf a))) 1e-15)
    (check-= (Fresnel-C a) (bigfloat->flonum (bfFresnel-C (bf a))) 1e-15))

  (test 0.1)
  (test 0.5)
  (test 1.0)
  (test 1.5)
  (test 1.7)
  (test 2.2)
  (test 2.5)
  (test 3.3)
  (test 4.4)
  ;(test 15.15)
  (bf-precision 1024)
  (test 15.15)

  (check-equal? (Fresnel-S -1)(-(Fresnel-S 1)))
  (check-equal? (Fresnel-C -1)(-(Fresnel-C 1)))
  (check-equal? (Fresnel-S -5)(-(Fresnel-S 5)))
  (check-equal? (Fresnel-C -5)(-(Fresnel-C 5)))
    
  )
