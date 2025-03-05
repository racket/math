#lang racket/base

(require racket/flonum "../../base.rkt")

(provide (all-defined-out))

(define (eqv x)
  (cond [(flonum? x) fl=]
        #;[(real? x) =]
        #;[(float-complex? x) =]
        [else =]))
(define (sub x)
  (cond [(flonum? x) fl-]
        #;[(real? x) -]
        #;[(float-complex? x) -]
        [else -]))
(define (div x)
  (cond [(flonum? x) fl/]
        #;[(real? x) /]
        #;[(float-complex? x) /]
        [else /]))
(define add
  (let ([fc+ (λ x* (apply + 0.0 x*))])
    (λ (x)
      (cond [(flonum? x) fl+]
            #;[(real? x) +]
            [(float-complex? x) fc+]
            [else +]))))
(define mul
  (let ([fc* (λ x* (apply * 1.0 x*))])
    (λ (x)
      (cond [(flonum? x) fl*]
            #;[(real? x) *]
            [(float-complex? x) fc*]
            [else *]))))
(define add*
  (let ([fc+ (λ x* (apply + 0.0+0.0i x*))])
    (λ (x)
      (cond [(flonum? x) fl+]
            #;[(real? x) +]
            [(float-complex? x) fc+]
            [else +]))))
(define mul*
  (let ([fc* (λ x* (apply * 1.0+0.0i x*))])
    (λ (x)
      (cond [(flonum? x) fl*]
            #;[(real? x) *]
            [(float-complex? x) fc*]
            [else *]))))
