#lang typed/racket/base

(require racket/promise
         math/distributions
         math/flonum (only-in racket/math nan?)
         typed/rackunit)

(define tol 1e-12)

;; size = 1 is the geometric distribution
(for ([k (in-range 0 20)])
  (check-= (pdf (negative-binomial-dist 1 0.3) k) (pdf (geometric-dist 0.3) k) tol)
  (check-= (cdf (negative-binomial-dist 1 0.3) k) (cdf (geometric-dist 0.3) k) tol))

;; hand-computed: Γ(5.5)/(3! Γ(2.5)) * 0.4^2.5 * 0.6^3
(define D (negative-binomial-dist 2.5 0.4))
(check-equal? (negative-binomial-dist-size D) 2.5)
(check-equal? (negative-binomial-dist-prob D) 0.4)
(check-= (pdf D 3) 0.14344091466523773 tol)
(check-= (cdf D 3) 0.5558019215511943 tol)
(check-= (exp (pdf D 3 #t)) (pdf D 3) tol)

;; integer size: P[NB(r,p) <= k] = P[Binomial(r+k,p) >= r]
(for* ([r (in-list '(1 3 7))] [k (in-range 0 15)])
  (check-= (cdf (negative-binomial-dist r 0.35) k)
           (cdf (binomial-dist (+ r k) 0.35) (- r 1) #f #t)
           tol))

;; non-integers and out-of-support values
(check-equal? (pdf D -1) 0.0)
(check-true (nan? (pdf D 1.5)))
(check-equal? (cdf D -1) 0.0)

;; inverse cdf round trip, and median
(for ([k (in-range 0 20)])
  (check-equal? (inv-cdf D (cdf D k)) (->fl k)))
(check-equal? (force (ordered-dist-median D)) (inv-cdf D 1/2))

;; undefined distributions
(check-true (nan? (pdf (negative-binomial-dist -1 0.5) 2)))
(check-true (nan? (pdf (negative-binomial-dist 2 1.5) 2)))

;; samples are nonnegative integers
(check-true (andmap (λ: ([x : Real]) (and (>= x 0) (integer? x))) (sample D 50)))