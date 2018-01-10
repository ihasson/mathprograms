(require xrepl)

(define (double x)
  (if (real? x)
    (real->double-flonum x)
    x))

(define (bisec-method a b ε n f)
  (let ((p (/ (+ a b) 2)))
    (if (and (> (abs (f p)) ε) (< n 1000))
      (if (> (* (f a) (f p)) 0)
        (bisec-method p b ε (+ n 1) f)
        (bisec-method a p ε (+ n 1) f))
      (begin 
        (displayln 
          (list "p =" (double p) " , n =" n " , f(p) =" (double(f p)) ))
        p))))

(define g (lambda (x) (- (* x x) 1)))

(define f3 (lambda (x) (+ (expt x 3) (* -7 (expt x 2)) (* 14 x) -6)))

