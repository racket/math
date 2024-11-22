#lang typed/racket/base

(require racket/performance-hint
         racket/promise
         "../../flonum.rkt"
         "../unsafe.rkt"
         "impl/student-t.rkt"
         "utils.rkt")

(provide Student-T-Dist
         student-t-dist
         student-t-dist-freedom    ; ν
         student-t-dist-mean       ; μ
         student-t-dist-scale      ; σ
         )


;; ===================================================================================================
;; Distribution object

(define-real-dist: student-t-dist Student-T-Dist
  student-t-dist-struct ([freedom : Positive-Float] [mean : Real] [scale : Real])) ; ν, μ, σ

(begin-encourage-inline
  (: student-t-dist (case-> (Real           -> Student-T-Dist)
                            (Real Real Real -> Student-T-Dist)))
  (define (student-t-dist ν [μ 0] [σ 1])
    (let ([ν* (fl ν)] [μ* (fl μ)] [σ* (fl σ)])
      (cond
        [(<= ν* 0) (raise-argument-error 'student-t-dist (format "Positive Real (~a < ν)" +min.0) 0 ν*)]
        [else
         (define std? (and (= μ 0.) (= σ 1.)))
         (define-values (pdf cdf inv-cdf sample)
           (if std?
               (values (make-student-t-pdf ν)
                       (make-student-t-cdf ν)
                       (make-student-t-inverse-cdf ν)
                       (case-lambda:
                         [()               (unsafe-flvector-ref (flstudent-t-sample ν* 1) 0)]
                         [([n : Integer])  (flvector->list (flstudent-t-sample ν* n))]))
               (values (make-student-t-pdf ν μ σ)
                       (make-student-t-cdf ν μ σ)
                       (make-student-t-inverse-cdf ν μ σ)
                       (case-lambda:
                         [()               (unsafe-flvector-ref (flstudent-t-sample ν* μ* σ* 1) 0)]
                         [([n : Integer])  (flvector->list (flstudent-t-sample ν* μ* σ* n))]))))
         (define median (delay μ*))
         ((inst student-t-dist-struct Real Flonum) pdf sample cdf inv-cdf
                                                   -inf.0 +inf.0
                                                   median ν* μ σ)]))))
