#lang racket/base

(require typed/untyped-utils)

(require (except-in "private/functions/gamma.rkt" gamma)
         (except-in "private/functions/log-gamma.rkt" log-gamma)
         (except-in "private/functions/beta.rkt" beta log-beta)
         (except-in "private/functions/erf.rkt" erf erfc)
         (except-in "private/functions/lambert.rkt" lambert)
         (except-in "private/functions/zeta.rkt" eta zeta)
         (except-in "private/functions/hurwitz-zeta.rkt" hurwitz-zeta)
         (except-in "private/functions/psi.rkt" psi0)
         "private/functions/incomplete-gamma.rkt"
         "private/functions/incomplete-beta.rkt"
         "private/functions/stirling-error.rkt"
         (except-in "private/functions/fresnel.rkt" Fresnel-S Fresnel-RS Fresnel-C Fresnel-RC))

(require/untyped-contract
 "private/functions/gamma.rkt"
 [gamma  (Number -> Number)])

(require/untyped-contract
 "private/functions/log-gamma.rkt"
 [log-gamma  (Number -> Number)])

(require/untyped-contract
 "private/functions/psi.rkt"
 [psi0  (Number -> Number)])

(require/untyped-contract
 "private/functions/beta.rkt"
 [beta  (Real Real -> Real)]
 [log-beta  (Real Real -> Real)])

(require/untyped-contract
 "private/functions/erf.rkt"
 [erf   (Number -> Number)]
 [erfc  (Real -> Real)])
(require/untyped-contract
 "private/functions/fresnel.rkt"
 [Fresnel-S  (Number -> Number)]
 [Fresnel-RS (Number -> Number)]
 [Fresnel-C  (Number -> Number)]
 [Fresnel-RC (Number -> Number)])


(require/untyped-contract
 "private/functions/lambert.rkt"
 [lambert  (Real -> Real)])

(require/untyped-contract
 "private/functions/zeta.rkt"
 [eta   (Real -> Real)]
 [zeta  (Real -> Real)])

(require/untyped-contract
 "private/functions/hurwitz-zeta.rkt"
 [hurwitz-zeta  (Real Real -> Real)])

(provide (all-from-out
          "private/functions/gamma.rkt"
          "private/functions/log-gamma.rkt"
          "private/functions/beta.rkt"
          "private/functions/erf.rkt"
          "private/functions/lambert.rkt"
          "private/functions/zeta.rkt"
          "private/functions/hurwitz-zeta.rkt"
          "private/functions/psi.rkt"
          "private/functions/incomplete-gamma.rkt"
          "private/functions/incomplete-beta.rkt"
          "private/functions/stirling-error.rkt"
          "private/functions/fresnel.rkt")
         gamma
         log-gamma
         psi0
         beta log-beta
         erf erfc
         lambert
         eta zeta
         hurwitz-zeta
         Fresnel-S Fresnel-RS
         Fresnel-C Fresnel-RC)
