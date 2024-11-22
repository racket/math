#lang typed/racket/base

(require (for-syntax typed/racket/base)
         math/distributions
         math/flonum
         typed/rackunit)

(define-syntax (check~= stx)
  (syntax-case stx ()
    [(_ v1 v2 i) (syntax/loc stx (check-= v1 v2 (* (abs v1) i epsilon.0)))]))

(define T (student-t-dist 2))
(check~= (pdf T 1)  (/ 1 (* 3 (sqrt 3)))            1)
(check~= (cdf T 1)  (+ 1/2 (/ 1 (* 2 (sqrt 3))))    1)
(check~= (inv-cdf T (+ 1/2 (/ 1 (* 2 (sqrt 3))))) 1 1)

(define pi (angle -1))
; N[PDF[StudentTDistribution[1], 0], 30]
(check~= (pdf (student-t-dist 1) 0)      0.3183098861837907       1)
(check~= (pdf (student-t-dist 2) 1)      (/ 1 (* 3 (sqrt 3)))     1)
(check~= (pdf (student-t-dist 3) 2)      (/ (* 6 (sqrt 3)) 49 pi) 1)
; generalized
;; x=1 | μ=2 | σ=1 => y=(x-μ)/σ = -1 => == (pdf (student-t-dist 2) ±1)
(check~= (pdf (student-t-dist 2 2 1) 1)      (/ 1 (* 3 (sqrt 3))) 1)
(check~= (pdf (student-t-dist 1 2 2) 0)      0.07957747154594767  1)
(check~= (pdf (student-t-dist 1 2 2) 1)      0.12732395447351627  1)
(check~= (pdf (student-t-dist 5 3 4) 0)      0.06892452901798418  1)
(check~= (pdf (student-t-dist 5 3 4) 1)      0.08197963283068663  1)
; Log space
(check~= (pdf (student-t-dist 2) 1 #t)   (fllog (pdf (student-t-dist 2) 1)) 1)

; N[CDF[StudentTDistribution[1], 0], 30]
(check-equal? (cdf (student-t-dist 1)  0)      0.5)
(check-equal? (cdf (student-t-dist 1)  1)      0.75)
(check~=      (cdf (student-t-dist 1) -1)      0.25                2)
(check-equal? (cdf (student-t-dist 2)  0)      0.5)
(check-equal? (cdf (student-t-dist 2)  1)      0.7886751345948129)
(check~=  (cdf (student-t-dist 2) -1)          0.2113248654051871  2)
; generalized
(check~= (cdf (student-t-dist 3 1 2) -1)       0.19550110947788532 1)
(check~= (cdf (student-t-dist 3 1 2)  0)       0.3257239824240755  1)
(check~= (cdf (student-t-dist 3 1 2)  1)       0.5                 1)
; Log space
(check~= (cdf (student-t-dist 2) 1 #t)   (fllog (cdf (student-t-dist 2) 1))  2)

;   N[InverseCDF[StudentTDistribution[2], 1/10], 30]
(define 2epsilon.0 (* 2. epsilon.0))
; Special case ν=1
(check-equal?  (inv-cdf (student-t-dist 1)  0)   -inf.0)
(check-equal?  (inv-cdf (student-t-dist 1)  1)   +inf.0)
(check-equal?  (inv-cdf (student-t-dist 1)  0.5) 0.) 
(check~=       (inv-cdf (student-t-dist 1)  0.1) -3.0776835371752536 1)
(check~=       (inv-cdf (student-t-dist 1)  0.9)  3.0776835371752536 1)
; Special case ν=2
(check-equal?  (inv-cdf (student-t-dist 2)  0)   -inf.0)
(check-equal?  (inv-cdf (student-t-dist 2)  1)   +inf.0)
(check-equal?  (inv-cdf (student-t-dist 2)  0.5) 0.) 
(check~=       (inv-cdf (student-t-dist 2)  0.1) -1.8856180831641267 1)
(check~=       (inv-cdf (student-t-dist 2)  0.9)  1.8856180831641267 1)
; Special case ν=4
(check-equal?  (inv-cdf (student-t-dist 4)  0)   -inf.0)
(check-equal?  (inv-cdf (student-t-dist 4)  1)   +inf.0)
(check-equal?  (inv-cdf (student-t-dist 4)  0.5) 0.)
(check-equal?  (inv-cdf (student-t-dist 4.)  0.5) 0.)
(check~=       (inv-cdf (student-t-dist 4)  0.1) -1.5332062740589438 1)
(check~=       (inv-cdf (student-t-dist 4)  0.9)  1.5332062740589438 2)
; General case
(check-equal?  (inv-cdf (student-t-dist 3)  0.5) 0.) 
(check~=       (inv-cdf (student-t-dist 3)  0.1) -1.6377443536962102 1)
(check~=       (inv-cdf (student-t-dist 3)  0.9)  1.6377443536962102 2)
(check-equal?  (inv-cdf (student-t-dist 5)  0)   -inf.0)
(check-equal?  (inv-cdf (student-t-dist 5)  1)   +inf.0)
(check-equal?  (inv-cdf (student-t-dist 5)  0.5) 0.) 
(check~=       (inv-cdf (student-t-dist 5)  0.1) -1.475884048824481 2)
(check~=       (inv-cdf (student-t-dist 5)  0.9)  1.475884048824481 3)
; Three parameters
(check~=       (inv-cdf (student-t-dist 3 1 2)  0.5) 1                   1)
(check~=       (inv-cdf (student-t-dist 3 1 2)  0.1) -2.2754887073924204 1)
(check~=       (inv-cdf (student-t-dist 3 1 2)  0.9)  4.27548870739242   1)
; Log space
(check~= (inv-cdf (student-t-dist 2 ) (log .25) #t)   (inv-cdf (student-t-dist 2 ) .25)  1)
(check~= (inv-cdf (student-t-dist 2.) (log .25) #t)   (inv-cdf (student-t-dist 2.) .25)  1)

