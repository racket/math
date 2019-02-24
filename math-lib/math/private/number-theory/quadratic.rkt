#lang typed/racket/base

(require racket/list
         "types.rkt")

(provide quadratic-complex-solutions
         quadratic-solutions
         quadratic-integer-solutions
         quadratic-natural-solutions)

(: quadratic-complex-solutions  : Complex Complex Complex -> (Listof Complex))
(define (quadratic-complex-solutions a b c)
  ; return list of solutions to a a x^2 + b x + c = 0
  (let ([d (- (* b b) (* 4 a c))])
    (if (= d 0)
        (make-list 2 (/ b (* -2 a)))
        (let* ([sqrt-d (sqrt d)]
               [sign (if (>= 0 (real-part (* (conjugate b) sqrt-d))) 1 -1)]
               [q (/ (+ b (* sign sqrt-d)) -2)])
          (list (/ q a) (/ c q))))))

(: quadratic-solutions : Complex Complex Complex -> (Listof Real))
(define (quadratic-solutions a b c)
  (filter real? (quadratic-solutions a b c)))

(: quadratic-integer-solutions : Complex Complex Complex -> (Listof Integer))
(define (quadratic-integer-solutions a b c)
  ; return list of integer solutions to a x^2 + b x + c = 0
  (filter exact-integer? (quadratic-solutions a b c)))

(: quadratic-natural-solutions : Complex Complex Complex -> (Listof Natural))
(define (quadratic-natural-solutions a b c)
  ; return list of nonnegative-integer solutions to a x^2 + b x + c = 0
  (filter natural? (quadratic-solutions a b c)))
