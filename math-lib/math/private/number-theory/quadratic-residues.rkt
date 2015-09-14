#lang typed/racket/base

(require racket/list
         "divisibility.rkt"
         "modular-arithmetic.rkt"
         "number-theory.rkt")

(provide quadratic-character quadratic-residue? jacobi-symbol)

; DEFINITION (Quadratic residue)
;   a in Un is a quadratic residue,
;   if there exists an s such that a=s^2 (mod n)
;   The number s is called a squre root of a modulo n.

; p is prime
(: quadratic-character : Integer Integer -> (U -1 0 1))
(define (quadratic-character a p)
  (cond [(a . < . 0)  (raise-argument-error 'quadratic-character "Natural" 0 a p)]
        [(p . <= . 0)  (raise-argument-error 'quadratic-character "Positive-Integer" 1 a p)]
        [else  (let ([l  (modular-expt a (quotient (- p 1) 2) p)])
                 (cond [(or (eqv? l 0) (eqv? l 1))  l]
                       [else  -1]))]))

(: quadratic-residue? : Integer Integer -> Boolean)
(define (quadratic-residue? a n)
  (cond [(a . < . 0)  (raise-argument-error 'quadratic-residue? "Natural" 0 a n)]
        [(n . <= . 0)  (raise-argument-error 'quadratic-residue? "Positive-Integer" 1 a n)]
        [else
         (let* ([ps     (prime-divisors n)]
                [odd-ps (if (= (first ps) 2)
                            (rest ps)
                            ps)])
           (and (andmap (λ: ([p : Natural]) 
                          (= (quadratic-character a p) 1))
                        odd-ps)
                (cond
                  [(divides? 8 n)  (= (modulo a 8) 1)]
                  [(divides? 4 n)  (= (modulo a 4) 1)]
                  [else            #t])))]))

(: jacobi-symbol : Integer Positive-Integer -> (U -1 0 1))
(define (jacobi-symbol a n)
  (unless (odd? n)
    (raise-argument-error 'jacobi "odd?" 1 a n))
  (cond
    [(= n 1) 1]
    [else
     (define prime-factors (factorize n))
     (let next ([factor (first prime-factors)] [remaining-factors (rest prime-factors)])
       (define qcap (quadratic-character a (first factor)))
       (cast (if (null? remaining-factors)
           (expt qcap (second factor))
           (* (expt qcap (second factor)) (next (first remaining-factors) (rest remaining-factors)))) (U -1 0 1)))]))