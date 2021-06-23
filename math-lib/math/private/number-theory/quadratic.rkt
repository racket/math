#lang typed/racket/base

(require racket/list
         (only-in racket/math conjugate)
         "types.rkt")

(provide complex-quadratic-solutions  ; complex coefficients
         quadratic-solutions          ; real coefficients
         quadratic-integer-solutions  ; coefficients
         quadratic-natural-solutions) ;coefficients


(: complex-quadratic-solutions  : Complex Complex Complex -> (Listof Complex))
(define (complex-quadratic-solutions a b c)
  ; Return list of solutions to a a x^2 + b x + c = 0
  ; where a,b and c are complex numbers.
  (let* ([d         (- (* b b) (* 4 a c))]
         [sqrt-d    (sqrt d)]
         [-b-sqrt-d (- (- b) sqrt-d)]
         [2a        (* 2 a)])
    (cond
      [(= d 0) (list (/ b (- 2a)))]
      [(= c 0) (list (/ (- (- b) sqrt-d) 2a)
                     (/ (+ (- b) sqrt-d) 2a))]
      ; use the standard formula, unless -b and sqrt are almost equal
      [#t ; (> (magnitude -b-sqrt-d) 0.001)
       (list (/ -b-sqrt-d        2a)
             (/ (+ (- b) sqrt-d) 2a))]
      [else
       ; Note: Disabled for now.
       ;       There are cases where only one root needs to
       ;       use Muller's formula. But which one is it?
       ; Muller's formula:
       ;   x = 2c / ( -b +- sqrt(d) )
       (let* ([sign (if (>= 0 (real-part (* (conjugate b) sqrt-d)))
                        1 -1)]
              [q    (/ (+ b (* sign sqrt-d)) -2)])
         (list (/ q a) (/ c q)))])))

(: quadratic-solutions : Real Real Real -> (Listof Real))
(define (quadratic-solutions a b c)
  ; return list of solutions to a a x^2 + b x + c = 0
  (let ([d (- (* b b) (* 4 a c))])
    (cond
      [(< d 0) '()]
      [(= d 0) (list (/ b (* -2 a)))]
      [else
       (let ([sqrt-d (sqrt d)])
         (if (>= b 0)
             (list (/ (- (- b) sqrt-d) (* 2 a))
                   (/ (* 2 c) (- (- b) sqrt-d)))
             (list (/ (+ (- b) sqrt-d) (* 2 a))
                   (/ (* 2 c) (+ (- b) sqrt-d)))))])))

(: quadratic-integer-solutions : Real Real Real -> (Listof Integer))
(define (quadratic-integer-solutions a b c)
  ; return list of integer solutions to a x^2 + b x + c = 0
  (filter exact-integer? (quadratic-solutions a b c)))

(: quadratic-natural-solutions : Real Real Real -> (Listof Natural))
(define (quadratic-natural-solutions a b c)
  ; return list of nonnegative-integer solutions to a x^2 + b x + c = 0
  (filter natural? (quadratic-solutions a b c)))
