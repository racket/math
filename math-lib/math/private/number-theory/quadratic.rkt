#lang typed/racket/base

(require racket/list
         (only-in racket/math conjugate sgn nan?)
         "types.rkt"
         "../flonum/flonum-functions.rkt")

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
  (define ac-sqrt (* (flsqrt (real->double-flonum a))
                     (flsqrt (real->double-flonum c))))
  (define b/2 (/ b 2))

  (define sqrt-d
    (cond
     [(and (exact? a) (exact? b) (exact? c))
      ; If a/b/c are exact, we want to keep as much of the computation
      ; as possible exact, so the sqrt goes at the end
      (define sqrt-d-number (sqrt (- (* b/2 b/2) (* a c))))
      (if (real? sqrt-d-number)
          sqrt-d-number
          +nan.0)]
     [(= (sgn a) (sgn c))
      ; In this case we know that ac is positive so ac-sqrt has the
      ; right sign. In this case we use difference of squares, which
      ; allows us to do the sqrt operations without overflowing.
      ;
      ; So the plan is sqrt(|b/2| + ac-sqrt) sqrt(|b/2| - ac-sqrt)
      (let* ([term1 (- (abs b/2) ac-sqrt)]
             [term2 (+ (abs b/2) ac-sqrt)])
             ; But there is a danger here that ac-sqrt rounded itself, either
             ; up or down, and so the second term evaluates to zero, when it
             ; should not. In this case, ac-sqrt has an error of ac-sqrt *
             ; epsilon, and for the subtraction to yield zero we must have
             ; |b/2| ~~ ac-sqrt, so the overall magnitude of the error is b
             ; sqrt(epsilon). That's not good enough, so in this case we
             ; expand the series out by one tick. That brings the error to a
             ; small number of ulps, which is good enough.
        (* (flsqrt
            (if (= term1 0)
                ; In this case we just need the error term of ac-sqrt. We solve:
                ;
                ;   (ac-sqrt + err)^2 = a c
                ;   ac-sqrt^2 + 2 ac-sqrt err + err^2 = a c
                ;   err = (a c - ac-sqrt^2) / 2 ac-sqrt
                ;       = a / 2 * c / ac-sqrt - ac-sqrt/2
                ;
                ; In this derivation I'm dropping the err^2 term
                ; because it's too small to matter. The key is to
                ; avoid overflow in the final formula. The most
                ; important thing to remember here is that for c's and
                ; ac-sqrt's exponents to have the opposite signs (and
                ; thus for the final formula to result in overflow) we
                ; must have a's and c's exponents to have opposite
                ; signs, so a's and ac-sqrt's exponents must have the
                ; same sign.
                (if (or (and (> (abs a) 1.0) (> ac-sqrt 1.0))
                        (and (< (abs a) 1.0) (< ac-sqrt 1.0)))
                    (/ (- ac-sqrt (* (/ a ac-sqrt) c)) 2)
                    (/ (- ac-sqrt (* a (/ c ac-sqrt))) 2))
                term1))
           ; This one has no worries about cancellation because it's
           ; the sum of two positive values
           (flsqrt term2)))]
     [else
      ; In this case hypot is perfect.
      (flhypot (real->double-flonum b/2) ac-sqrt)]))

  ; return list of solutions to a a x^2 + b x + c = 0
  (cond
   [(nan? sqrt-d) 
    ; Use the standard a/c swap trick to avoid cancellation
    (if (< b 0)
        (list (/ (+ (- b/2) sqrt-d) a) (/ c (+ (- b/2) sqrt-d)))
        (list (/ (- (- b/2) sqrt-d) a) (/ c (- (- b/2) sqrt-d))))]
   [(zero? sqrt-d)
    (list (- (/ b/2 a)))]
   [else '()]))

(: quadratic-integer-solutions : Real Real Real -> (Listof Integer))
(define (quadratic-integer-solutions a b c)
  ; return list of integer solutions to a x^2 + b x + c = 0
  (filter exact-integer? (quadratic-solutions a b c)))

(: quadratic-natural-solutions : Real Real Real -> (Listof Natural))
(define (quadratic-natural-solutions a b c)
  ; return list of nonnegative-integer solutions to a x^2 + b x + c = 0
  (filter natural? (quadratic-solutions a b c)))
