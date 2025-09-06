#lang typed/racket/base

(require (for-syntax typed/racket/base)
         racket/promise
         math/distributions
         math/flonum
         typed/rackunit)

(define-syntax (check~= stx)
  (syntax-case stx ()
    [(_ v1 v2 i)
     (quasisyntax (let ([V1 v1][V2 v2])
                    (unless (eqv? V1 V2)
                      #,(syntax/loc stx (check-= V1 V2 (* (abs V1) i epsilon.0))))))]))

(define W1 (weibull-dist 2))
(check-equal? (weibull-dist-shape W1) 2.)
(check-equal? (weibull-dist-location W1) 0)
(check-equal? (weibull-dist-scale W1) 1)
(check~= (force (ordered-dist-median W1))  (flsqrt (fllog 2.))  1)

(check~= (pdf W1  1)  (* 2 (exp -1))  1)
(check~= (pdf W1 -1)  0  1)
(check~= (cdf W1  1)  (- 1 (exp -1))  1)
(check~= (inv-cdf W1 (- 1 (exp -1)))  1  1)
(check~= (inv-cdf W1 -1)  +nan.0  1)
(check~= (inv-cdf W1 2)  +nan.0  1)


(define W3 (weibull-dist 2 1/3 3))
(check-equal? (weibull-dist-shape W3) 2.)
(check-equal? (weibull-dist-location W3) 1/3)
(check-equal? (weibull-dist-scale W3) 3)
(check~= (force (ordered-dist-median W3))  (+ 1/3 (* 3 (flsqrt (fllog 2.))))  1)
(check~= (force (ordered-dist-median W3))  (inv-cdf W3 1/2)  1)

;; X = (/ (- x d) s)
(check~= (pdf W3  10/3) (* 2/3 (exp -1)) 1)
(check~= (pdf W3   0.3)    0             1) ; < 1/3 = 0.
(check~= (cdf W3  10/3) (- 1 (exp -1))   1)
(check~= (inv-cdf W3 (- 1 (exp -1))) (+ 1/3 (* 3 1))  1)
(check~= (inv-cdf W3  -1)  +nan.0  1)
(check~= (inv-cdf W3   2)  +nan.0  1)

(check-true (<= 1/3 (sample W3)))
(check-true (<= 1/3 (apply min (sample W3 20))))
(check-true (<= (apply max (sample (weibull-dist 1 1/3 -1) 20)) 1/3))
