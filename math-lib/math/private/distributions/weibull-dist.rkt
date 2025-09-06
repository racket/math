#lang typed/racket/base

(require racket/performance-hint
         racket/promise
         (only-in racket/math sgn)
         "../../flonum.rkt"
         "../unsafe.rkt"
         "impl/weibull.rkt"
         "utils.rkt")

(provide Weibull-Dist
         weibull-dist
         weibull-dist-shape    ; k
         weibull-dist-location ; d
         weibull-dist-scale    ; s (λ)
         )


;; ===================================================================================================
;; Distribution object

(define-real-dist: weibull-dist Weibull-Dist
  weibull-dist-struct ([shape : Positive-Float] [location : Real] [scale : Real])) ; k, d, s

(begin-encourage-inline
  (: weibull-dist (case-> (Real           -> Weibull-Dist)
                          (Real Real Real -> Weibull-Dist)))
  (define (weibull-dist k [d 0] [s 1])
    (define k* (fl k))
    (cond
      [(<= k* 0) (raise-argument-error 'weibull-dist (format "Positive Real (~a < ν)" +min.0) 0 k*)]
      [else
       (define std? (and (= s 0.) (= d 1.)))
       (define-values (pdf cdf inv-cdf sample)
         (if std?
             (let ([sampler (make-weibull-sampler k)])
               (values (make-weibull-pdf k)
                       (make-weibull-cdf k)
                       (make-weibull-inverse-cdf k)
                       (case-lambda:
                         [()               (unsafe-flvector-ref (sampler 1) 0)]
                         [([n : Integer])  (flvector->list (sampler n))])))
             (let ([sampler (make-weibull-sampler k d s)])
               (values (make-weibull-pdf k d s)
                       (make-weibull-cdf k d s)
                       (make-weibull-inverse-cdf k d s)
                       (case-lambda:
                         [()               (unsafe-flvector-ref (sampler 1) 0)]
                         [([n : Integer])  (flvector->list (sampler n))])))))
       (define median (delay (fl (+ d (* s (flexpt (fllog 2.) (/ k*)))))))
       ((inst weibull-dist-struct Real Flonum) pdf sample cdf inv-cdf
                                               (case (sgn s)
                                                 [(-1 0) -inf.0]
                                                 [else      0.0])
                                               (case (sgn s)
                                                 [(-1)    0.0]
                                                 [else +inf.0])
                                               median k* d s)])))
