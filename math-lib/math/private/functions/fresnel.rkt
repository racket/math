#lang typed/racket/base

#|
(natural)Fresnel Integrals:
 S(x)=∫sin(π/2t²) |0->x
 C(x)=∫cos(π/2t²) |0->x
or regular
 S(x)=∫sin(t²) |0->x
 C(x)=∫cos(t²) |0->x

floating point implementation:
  adaptation of the algorithm in the fresnl.c file from Cephes,
  but with the limit for the rational powerseries lowered to x<1.5625 (coming from 2.5625)
  above x>1.6 the error becomes too big (~1e-8 x~2.5)
Bigfloat implementation:
  adaptation of powerseries as shown on wikipedia

would like to extend this together with erf to the complex plane

|#

(require "../../base.rkt"
         "../../flonum.rkt"
         "erf.rkt")

(provide flFresnel-S complex-Fresnel-S Fresnel-S Fresnel-RS
         flFresnel-C complex-Fresnel-C Fresnel-C Fresnel-RC)

;------------------------------
;polinomials for the rational assymptotic aproximation for S & C
(define fn/d (make-quotient-flpolyfun
              ( 3.76329711269987889006E-20
                1.34283276233062758925E-16
                1.72010743268161828879E-13
                1.02304514164907233465E-10
                3.05568983790257605827E-8
                4.63613749287867322088E-6
                3.45017939782574027900E-4
                1.15220955073585758835E-2
                1.43407919780758885261E-1
                4.21543555043677546506E-1)
              ( 1.25443237090011264384E-20
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
(define gn/d (make-quotient-flpolyfun
              ( 1.86958710162783235106E-22
                8.36354435630677421531E-19
                1.37555460633261799868E-15
                1.08268041139020870318E-12
                4.45344415861750144738E-10
                9.82852443688422223854E-8
                1.15138826111884280931E-5
                6.84079380915393090172E-4
                1.87648584092575249293E-2
                1.97102833525523411709E-1
                5.04442073643383265887E-1)
              ( 1.86958710162783236342E-22
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
(define sn/d
  (make-quotient-flpolyfun ( 3.18016297876567817986E11
                            -4.42979518059697779103E10
                             2.54890880573376359104E9
                            -6.29741486205862506537E7
                             7.08840045257738576863E5
                            -2.99181919401019853726E3
                             0.0)
                           ( 6.07366389490084639049E11
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
     (fl* (fl* X2 x) (sn/d X4))]
    [else
     (define t (fl* pi (fl* x x)))
     (define t/2 (fl/ t 2.0))
     (define U (fl/ 1.0 (fl* t t)))
     (define f (fl- 1.0 (fl* U (fn/d U))))
     (define g (fl/ (gn/d U) t))
     (fl- 0.5 (fl/ (fl+ (fl* f (flcos t/2))(fl* g (flsin t/2)))
                   (fl* pi x)))]))

(define (complex-Fresnel-S [z : Number]) : Number
  (define z* (* z (sqrt (/ pi 4))))
  (* (/ 1+i 4)
     (-       (complex-erf (* z* 1+i))
        (* +i (complex-erf (* z* 1-i))))))

;------------------------------
;polinomials for the rational powerseries aproximation for C
(define cn/d
  (make-quotient-flpolyfun ( 9.99999999999999998822E-1
                            -2.05525900955013891793E-1
                             1.88843319396703850064E-2
                            -6.45191435683965050962E-4
                             9.50428062829859605134E-6
                            -4.98843114573573548651E-8
                             0.0)
                           ( 1.00000000000000000118E0
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
     (fl* x (cn/d X4))]
    [else
     (define t (fl* pi (fl* x x)))
     (define t/2 (fl/ t 2.0))
     (define U (fl/ 1.0 (fl* t t)))
     (define f (fl- 1.0 (fl* U (fn/d U))))
     (define g (fl/ (gn/d U) t))
     (fl+ 0.5 (fl/ (fl- (fl* f (flsin t/2))(fl* g (flcos t/2)))
                   (fl* pi x)))]))

(define (complex-Fresnel-C [z : Number]) : Number
  (define z* (* z (sqrt (/ pi 4))))
  (* (/ 1-i 4)
     (+       (complex-erf (* z* 1+i))
        (* +i (complex-erf (* z* 1-i))))))

;------------------------------
(: Fresnel-S (case-> (Zero -> Zero)
                     (Flonum -> Flonum)
                     (Real -> (U Zero Flonum))
                     (Number -> Number)))
(define (Fresnel-S z)
  (define x (if (= (imag-part z) 0) (real-part z) z))
  (cond
    [(eqv? z 0) 0]
    [(flonum? x)(flFresnel-S x)]
    [(real? x)(flFresnel-S (fl x))]
    [else (complex-Fresnel-S z)]))
(: Fresnel-C (case-> (Zero -> Zero)
                     (Flonum -> Flonum)
                     (Real -> (U Zero Flonum))
                     (Number -> Number)))
(define (Fresnel-C z)
  (define x (if (= (imag-part z) 0) (real-part z) z))
  (cond
    [(eqv? z 0) 0]
    [(flonum? x)(flFresnel-C x)]
    [(real? x)(flFresnel-C (fl x))]
    [else (complex-Fresnel-C z)]))
(: Fresnel-RS (case-> (Zero -> Zero)
                      (Flonum -> Flonum)
                      (Real -> (U Zero Flonum))
                      (Number -> Number)))
(define (Fresnel-RS z)
  (define x (if (= (imag-part z) 0) (real-part z) z))
  (cond
    [(eqv? z 0) 0]
    [(flonum? x)(* (flsqrt (fl/ pi 2.0))(flFresnel-S (* (flsqrt (fl/ 2.0 pi)) x)))]
    [(real? x)  (* (flsqrt (fl/ pi 2.0))(flFresnel-S (* (flsqrt (fl/ 2.0 pi)) (fl x))))]
    [else (* (sqrt (/ pi 2))(complex-Fresnel-S (* (sqrt (/ 2 pi)) z)))]))
(: Fresnel-RC (case-> (Zero -> Zero)
                      (Flonum -> Flonum)
                      (Real -> (U Zero Flonum))
                      (Number -> Number)))
(define (Fresnel-RC z)
  (define x (if (= (imag-part z) 0) (real-part z) z))
  (cond
    [(eqv? z 0) 0]
    [(flonum? x)(* (flsqrt (fl/ pi 2.0))(flFresnel-C (* (flsqrt (fl/ 2.0 pi)) x)))]
    [(real? x)  (* (flsqrt (fl/ pi 2.0))(flFresnel-C (* (flsqrt (fl/ 2.0 pi)) (fl x))))]
    [else (* (sqrt (/ pi 2))(complex-Fresnel-C (* (sqrt (/ 2 pi)) z)))]))

;------------------------------
(module* bfFresnel #f
  (require math/bigfloat)
  (provide bfFresnel-S bfFresnel-RS bfFresnel-C bfFresnel-RC)
  ;this implementation is potentially exact
  ;but for large x (> 5!) it is really slow and needs a lot of bits in bf-precision (~2x²)!
  (define (precision-check [a : Bigfloat][maxp : Integer]) : Integer
    (define p (bf-precision))
    (define a2 (expt (bigfloat->real a) 2))
    (define min-precision-needed (* 2 a2))
    (define expected-precision-loss (* 4/5 a2))
    (define mp (+ 5 (round (inexact->exact (max min-precision-needed (+ p expected-precision-loss))))))
    (cond
      [(<= mp maxp) mp]
      [else
       (error (format "bfFresnel: calculation aborted
 Minimum precision needed for calculating ~a... is ~a
 This is more than the maximum allowed calculating precision ~a->~a."
                      (bigfloat->flonum a) mp p maxp))]))
  (define (bfFresnel-RS [x : Bigfloat][maxp : Integer (* 2 (bf-precision))])
    (define p (bf-precision))
    (bfcopy
     (parameterize ([bf-precision (precision-check x maxp)])
       (define X3 (bfexpt x (bf 3)))
       (define X4 (bfexpt x (bf 4)))
       (define prsn (bfexpt (bf 1/2) (bf p)))
       (define-values (s l)
         (for/fold : (Values Bigfloat Bigfloat)
           ([s : Bigfloat (bf/ X3 (bf 3))]
            [l : Bigfloat (bf/ X3 (bf 3))])
           ([n (in-naturals 1)]
            #:break (bf< (bfabs (bf/ l s)) prsn))
           (define l+ (bf* (bf -1) X4 (bf* l (bf (/ (- (* 4 n) 1)
                                                    (* 2 n)(+ (* 2 n) 1)(+ (* 4 n) 3))))))
           (values (bf+ s l+) l+)))
       s)))
  (define (bfFresnel-S [x : Bigfloat][maxp : Integer (* 2 (bf-precision))])
    (bf* (bfsqrt (bf/ (bf 2) pi.bf))
         (bfFresnel-RS (bf* (bfsqrt (bf/ pi.bf (bf 2))) x) maxp)))
  
  (define (bfFresnel-RC [x : Bigfloat][maxp : Integer (* 2 (bf-precision))])
    (define p (bf-precision))
    (bfcopy
     (parameterize ([bf-precision (precision-check x maxp)])
       (define X4 (bfexpt x (bf 4)))
       (define prsn (bfexpt (bf 1/2) (bf (bf-precision))))
       (define-values (s l)
         (for/fold : (Values Bigfloat Bigfloat)
           ([s : Bigfloat x]
            [l : Bigfloat x])
           ([n (in-naturals 1)]
            #:break (bf< (bfabs (bf/ l s)) prsn))
           (define l+ (bf* (bf -1) X4 (bf* l (bf (/ (- (* 4 n) 3)
                                                    (* 2 n)(- (* 2 n) 1)(+ (* 4 n) 1))))))
           (values (bf+ s l+) l+)))
       s)))
  (define (bfFresnel-C [x : Bigfloat][maxp : Integer (* 2 (bf-precision))])
    (bf* (bfsqrt (bf/ (bf 2) pi.bf))
         (bfFresnel-RC (bf* (bfsqrt (bf/ pi.bf (bf 2))) x) maxp)))
    )

#;(module* test #f
    ;some functions to check the function.
    ;commented out, because (partly) put in function-tests from math-test
  (require math/bigfloat
           typed/rackunit)
  (require (submod "..")
           (submod ".." bfFresnel))
  
  (define (test [a : Flonum] #:ε [ε : Flonum 1e-15])
    (check-= (Fresnel-S  a) (bigfloat->flonum (bfFresnel-S  (bf a))) ε (format "(Fresnel-S ~a)" a))
    (check-= (Fresnel-RS a) (bigfloat->flonum (bfFresnel-RS (bf a))) ε (format "(Fresnel-RS ~a)" a))
    (check-= (Fresnel-C  a) (bigfloat->flonum (bfFresnel-C  (bf a))) ε (format "(Fresnel-C ~a)" a))
    (check-= (Fresnel-RC a) (bigfloat->flonum (bfFresnel-RC (bf a))) ε (format "(Fresnel-RC ~a)" a)))

  (test 0.1)
  (test 0.5)
  (test 1.0)
  (test 1.5)
  (test 1.7)
  (test 2.2)
  (test 2.5)
  (test 3.3)
  (test 4.4)
  ;(test 15.15);<-bf-precision to low (standard=128)
  (parameterize ([bf-precision 1024])
    (test 12.087951096163412);even though this is ok for C
    (test 15.15 #:ε 1e-14);this fails for RC @1e-15
    )
  (check-equal? (Fresnel-S -1)(-(Fresnel-S 1)))
  (check-equal? (Fresnel-C -1)(-(Fresnel-C 1)))
  (check-equal? (Fresnel-S -5)(-(Fresnel-S 5)))
  (check-equal? (Fresnel-C -5)(-(Fresnel-C 5)))

  (check-= (magnitude
            (/ (complex-Fresnel-S 1+i)
               -2.0618882191948404680807165366857086008159083237378680520+2.0618882191948404680807165366857086008159083237378680520i))
           1 1e-12)
  (check-= (magnitude
            (/ (complex-Fresnel-S 5+0.2i)
               0.47365635370953447150430003290670950910437504109567910708+0.73462461062246762668741695076291749315532498862606992695i))
           1 1e-12)
  (check-= (magnitude
            (/ (complex-Fresnel-S -8-25i)
               4.33491319138289340482459027250885195089165476877024e270+1.38537242126661103439768955547584123133806956947632e270i))
           1 1e-12)
  (check-= (magnitude
            (/ (complex-Fresnel-C 1+i)
               2.55579377810243902463452238835219584215662360420358429635+2.55579377810243902463452238835219584215662360420358429635i))
           1 1e-12)
  (check-= (magnitude
            (/ (complex-Fresnel-C 5+0.2i)
               1.237351377588089955209810162536250644901272375752411527+0.02602758966318992966794130566838646516498797340046208293i))
           1 1e-12)
  (check-= (magnitude
            (/ (complex-Fresnel-C -8-25i)
               1.38537242126661103439768955547584123133806956947632e270-4.33491319138289340482459027250885195089165476877024e270i))
           1 1e-12)

  #;(let ()
    (local-require plot)
    (define (mk)(* 100 (random)))
    (define L
      (for/list : (Listof (List Flonum Flonum))
        ([i (in-range 5000)])
        (define a (mk))
        (list a
              (- (Fresnel-S a)
                 (bigfloat->flonum
                  (bfFresnel-S (bf a) 50000))))))
    (define M+ (apply max ((inst map Real (List Flonum Flonum)) (inst cadr Real) L)))
    (define M- (apply min ((inst map Real (List Flonum Flonum)) (inst cadr Real) L)))
    (define K : (HashTable Real Real) (make-hash))
    (for ([i (in-list L)])
      (hash-update! K (if (= (cadr i) 0) -33 (inexact->exact (round (/ (log (abs (cadr i)))(log 10))))) add1 (λ () 0)))
    (plot (points L)
          #:y-max (* 1.2 M+)
          #:y-min (* 1.2 M-)
          #:width 1000)
    (plot (discrete-histogram
           ((inst sort (List Real Real))
            ((inst map (List Real Real)(Pairof Real Real))
             (λ ([x : (Pairof Real Real)]) (list (car x)(cdr x)))
             (hash->list K))
            < #:key car))))
    
  )