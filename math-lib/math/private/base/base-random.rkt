#lang typed/racket/base

(require racket/fixnum)

(provide random-bits random-natural random-integer)

;; Random bits are taken in blocks of this size:
(define block-bits : Nonnegative-Fixnum 29)
(define block-size : Natural (arithmetic-shift 1 block-bits))

(: random-bits (case-> (Integer -> Natural)
                       (Integer Pseudo-Random-Generator -> Natural)))
(define (random-bits bits [prng (current-pseudo-random-generator)])
  (cond [(bits . < . 0)  (raise-argument-error 'random-bits "Non-Negative-Integer" bits)]
        [(bits . = . 0) 0]
        [(not (fixnum? bits))  (raise-argument-error 'random-bits "reasonably sized integer" bits)]
        [(bits . <= . block-bits)  (random (ann (fxlshift 1 bits) Nonnegative-Fixnum) prng)]
        [(bits . <= . (fx* 2 block-bits))
	 ;; this case is not always in the fixnum range on 32-bit platforms
         (define rem-bits (fx- bits block-bits))
         (bitwise-ior (arithmetic-shift (random block-size prng) rem-bits)
		      (random (ann (fxlshift 1 rem-bits) Nonnegative-Fixnum) prng))]
        [else
         (define max-blocks (assert (quotient bits block-bits) index?))
         (define rem-bits (remainder bits block-bits))
         (let: loop : Natural ([blocks : Nonnegative-Fixnum  0]
                               [r : Natural  (random (fxlshift 1 rem-bits) prng)])
           (cond [(blocks . fx< . max-blocks)
                  (loop (fx+ blocks 1)
                        (bitwise-ior (arithmetic-shift r block-bits)
                                     (random block-size prng)))]
                 [else  r]))]))

(define random-max 4294967087)

(: random-natural (case-> (Integer -> Natural)
                          (Integer Pseudo-Random-Generator -> Natural)))
;; Returns a random integer in the interval [0..n)
(define (random-natural n [prng (current-pseudo-random-generator)])
  (cond
    [(n . <= . 0)  (raise-argument-error 'random-natural "Positive-Integer" n)]
    [(n . <= . random-max)  (random n prng)]
    [else
     (define bits (integer-length (- n 1)))
     (let loop ()
       (define r (random-bits bits prng))
       (if (r . >= . n) (loop) r))]))

(: random-integer (case-> (Integer Integer -> Integer)
                          (Integer Integer Pseudo-Random-Generator -> Integer)))
(define (random-integer a b [prng (current-pseudo-random-generator)])
  (let ([a  (min a b)] [b  (max a b)])
    (+ a (random-natural (- b a) prng))))
