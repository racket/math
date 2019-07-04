#lang typed/racket

(require math/base
         math/special-functions
         math/flonum
         typed/rackunit)

(: ε (Parameterof Float))
(define ε (make-parameter epsilon.0))

(: relative-error<= : Float Float Float -> Boolean)
(define (relative-error<= x correct epsilon)
  (<= (relative-error x correct) epsilon))

;; ---------------------------------------------------------------------------------------------------
;; expm1

(ε (* 2 epsilon.0))

(check-= (flexpm1 1.)      1.71828182845904523536028747135 (ε))
(check-= (flexpm1 1e-5)    0.000010000050000166667083 (ε))
(check-= (flexpm1 1e-10)   1.00000000005000000000166666667e-10 (ε))
(check-= (flexpm1 1e-15)   1.00000000000000050000000000000e-15 (ε))
(check-= (flexpm1 1e-15)   1.00000000000000050000000000000e-15 (ε))
(check-= (flexpm1 -1.)    -0.632120558828557678404476229839 (ε))
(check-= (flexpm1 -1e-10) -9.99999999950000000001666666667e-11 (ε))
(check-= (flexpm1 0.0)     0.0 (ε))

;; ---------------------------------------------------------------------------------------------------
;; gamma

(ε epsilon.0)

(check-equal? (flgamma -0.0) -inf.0)
(check-equal? (flgamma 0.0) +inf.0)
(check-equal? (flgamma 1.0) 1.0)
(check-equal? (flgamma 2.0) 1.0)
(check-equal? (flgamma 3.0) 2.0)
(check-equal? (flgamma 4.0) 6.0)
(check-= (flgamma -21.5) 1.31844491832155110297694106059e-20 (ε))
(check-true (relative-error<= (flgamma 1e-15) 9.9999999999999942278433509847e14 (ε)))
(check-true (relative-error<= (flgamma 142.5) 2.25990910998653224305124991671e244 (ε)))
(check-equal? (flgamma 172.0) 1.24101807021766782342484052410e309) ; = +inf.0
(check-true (relative-error<= (flgamma -1e-5) -100000.577225555552235029678062 (ε)))

(check-exn exn:fail:contract? (λ () (gamma 0)))
(check-equal? (gamma 4) 6)
(check-equal? (gamma 4.0) (flgamma 4.0))

;; ---------------------------------------------------------------------------------------------------
;; hypot

(ε (* 2 epsilon.0))
  
(check-equal? (flhypot 0. 0.) 0.)
(check-equal? (flhypot 3. 0.) 3.)
(check-equal? (flhypot 3. 4.) 5.)
(check-equal? (flhypot -3. -4.) 5.)
(check-equal? (flhypot 4. 3.) 5.)
(check-= (flhypot 2. 3.) 3.6055512754639892931 (ε))
(check-= (flhypot 3. 2.) 3.6055512754639892931 (ε))
  
;; ---------------------------------------------------------------------------------------------------
;; logs, exponents, bases, etc.

(check-equal? (fllog/base 2.0 (expt 2.0 5.0)) 5.0)

(check-false (power-of-two? 3))
(check-true (power-of-two? 2))
(check-true (power-of-two? 1))
(check-true (power-of-two? 1/2))
(check-false (power-of-two? 0))
(check-false (power-of-two? -1))

;; ---------------------------------------------------------------------------------------------------
;; Fresnel integrals
(ε (* 2 epsilon.0)); real maximum error (* 20 epsilon.0) ~ 5e-15

;check if function is odd
(check-equal? (Fresnel-S -1)  (- (Fresnel-S 1)))
(check-equal? (Fresnel-RS -1) (- (Fresnel-RS 1)))
(check-equal? (Fresnel-C -1)  (- (Fresnel-C 1)))
(check-equal? (Fresnel-RC -1) (- (Fresnel-RC 1)))

;check some values calculated with bfFresnel (submod "math/private/functions/fesnel.rkt" bfFresnel)
(check-= (Fresnel-S  1) 0.4382591473903547660767566966251526374943 (ε))
(check-= (Fresnel-S  2) 0.3434156783636982421953008159580684568872 (ε))
(check-= (Fresnel-S  5) 0.4991913819171168867519283804659916406959 (ε))
(check-= (Fresnel-S 10) 0.4681699785848822404033511108106786585115 (ε))
(check-= (Fresnel-S 20) 0.4840845359259538927147542448558347167707 (ε))
(check-= (Fresnel-S 50) 0.4936338025859387414532682397988025642247 (ε))

(check-= (Fresnel-C  1) 0.779893400376822829474206413652690136631  (ε))
(check-= (Fresnel-C  2) 0.4882534060753407545002235033572610376925 (ε))
(check-= (Fresnel-C  5) 0.5636311887040122311021074044130139830057 (ε))
(check-= (Fresnel-C 10) 0.4998986942055157236141518477355803595637 (ε))
(check-= (Fresnel-C 20) 0.4999873349723443881870062136976602164476 (ε))
(check-= (Fresnel-C 50) 0.4999991894307279679558101639817919070024 (ε))

