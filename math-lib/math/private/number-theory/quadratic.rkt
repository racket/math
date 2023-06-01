#lang typed/racket/base

(require racket/list
         (only-in racket/math conjugate sgn nan?)
         "types.rkt"
         "../flonum/flonum-functions.rkt"
         "../flonum/flonum-more-functions.rkt")

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

;;; This algorithm for complex quadratics is quite complex. A full
;;; description can be found here:
;;; https://pavpanchekha.com/blog/accurate-quadratic.html
(: quadratic-discriminant : Flonum Flonum Flonum -> Real)
(define (quadratic-discriminant a b c)
  "Compute sqrt(b^2 - a c) with very high accuracy"
  (let ([aa (abs a)] [ab (abs b)] [ac (abs c)]
        [sa (> a 0)] [sb (> b 0)] [sc (> c 0)])
    (define x (* (sqrt aa) (sqrt ac)))
    (cond
     ;; Otherwise we have two cases depending on the sign of a c
     ;; In this case a c is positive and we want sqrt(b^2 - a c)
     [(equal? sa sc)
      (define-values (ac/x ac/x.e)
        ;; Need to compute err(x) ~ (a b / x - x) / 2
        (if (equal? (> aa 1) (> x 1)) ; In this case do a / x first
            (let* ([a/x (/ aa x)]
                   [a/x.e (/ (flfma a/x x (- aa)) x)]
                   [ac/x (* a/x ac)])
              (values ac/x (+ (flfma a/x (- ac) ac/x) (* a/x.e ac))))
            (let* ([c/x (/ ac x)]
                   [c/x.e (/ (flfma c/x x (- ac)) x)]
                   [ac/x (* aa c/x)])
              (values ac/x (+ (flfma (- aa) c/x ac/x) (* c/x.e a))))))

      ;; Now we have d* = |b| - sqrt(ac)
      (define d* (- (- ab x) (/ (- (- ac/x x) ac/x.e) 2)))
      (cond
       [(> d* 0) (* (sqrt d*) (sqrt (+ ab x)))]
       [(= d* 0) 0]
       [else +nan.0])]
     ;; In this case, a c is negative and we want sqrt(b^2 + a c)
     [else
      (if (> ab x) ;; Standard sort trick for sqrt(x^2 + y^2)
          (let ([z (/ x ab)])
            (* ab (flsqrt (+ 1.0 (* z z)))))
          (let ([z (/ ab x)])
            (* x (flsqrt (+ (* z z) 1.0)))))])))

(: quadratic-solutions : Real Real Real -> (Listof Real))
(define (quadratic-solutions a b c)
  (define b/2 (/ b 2))
  (define sqrt-d
    (if (and (exact? a) (exact? b) (exact? c))
        (let ([d (- (* b/2 b/2) (* a c))])
          (if (< d 0) +nan.0 (sqrt d)))
        (quadratic-discriminant (fl a) (fl b/2) (fl c))))

  ; return list of solutions to a a x^2 + b x + c = 0
  (cond
   [(zero? a)
    (list (- (/ c b)))]
   [(nan? sqrt-d) 
    '()]
   [(zero? sqrt-d)
    (list (- (/ b/2 a)))]
   ; Use the standard a/c swap trick to avoid cancellation
   [(< b 0)
    (list (/ c (- sqrt-d b/2)) (/ (- sqrt-d b/2) a))]
   [else
    (list (/ (+ b/2 sqrt-d) (- a)) (/ (- c) (+ b/2 sqrt-d)))]))

(: quadratic-integer-solutions : Real Real Real -> (Listof Integer))
(define (quadratic-integer-solutions a b c)
  ; return list of integer solutions to a x^2 + b x + c = 0
  (filter exact-integer? (quadratic-solutions a b c)))

(: quadratic-natural-solutions : Real Real Real -> (Listof Natural))
(define (quadratic-natural-solutions a b c)
  ; return list of nonnegative-integer solutions to a x^2 + b x + c = 0
  (filter natural? (quadratic-solutions a b c)))
