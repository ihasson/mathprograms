(require math)

;(struct Sx-coeffs (a b c d))
(define S-table (make-hash))
(define (tridiag fn Xs coefficients-table)
  (let* 
    ([l₀ 1]
     [μ₀ 0]
     [z₀ 0]
     [x₀ (car Xs)]
     [x₁ (cadr Xs)]
     [h₀ (- x₁ x₀)]
     [a₀ (fn x₀)]
     [a₁ (fn x₁)]
     [next-c-l-z 
       (tridiagonal  
         fn (cddr Xs)
         x₁ x₀
         h₀
         μ₀ z₀
         a₁ a₀
         coefficients-table)]
     [c₁ (first next-c-l-z)]
     [c₀ (- z₀ (* μ₀ c₁))] 
     [b₀ (- (/ (- a₁ a₀) h₀) 
            (/ (* h₀ (+ c₁ (* 2 c₀))) 3))]
     [d₀ (/ (- c₁ c₀) (* 3 h₀))]
     [_ (hash-set! coefficients-table x₀ (list a₀ b₀ c₀ d₀))])
    coefficients-table))

(define tridiagonal 
  (λ (   fn Xs 
         curr-x prev-x 
         prev-h 
         prev-μ prev-z 
         curr-a prev-a
         coefficients-table)
    (if (empty? Xs) (list 0 1 0)
     ; (let ([l-n 1]
     ;       [z-n 0]
     ;       [c-n 0]);need to fill in 
      (let* 
        ([_ (displayln 
              (list Xs "curr-x" curr-x "prev-x" prev-x "prev-h" prev-h 
                    "prev-μ" prev-μ "prev-z" prev-z 
                    "curr-a" curr-a "prev-a" prev-a
                             coefficients-table))]
         [next-x (car Xs)]
         [next-a (fn next-x)]
         [curr-h (- next-x curr-x)]
         [curr-α (- (* (/ 3 curr-h) (- next-a prev-a)) 
                    (* (/ 3 prev-h) (- curr-a prev-a)))]
         [_ (displayln (list "currr-α" curr-α ))]
         [curr-l (- (* 2 (- next-x prev-x)) (* prev-h prev-μ))]
         [curr-μ (/ curr-h curr-l)]
         [_ (displayln (list "curr-l" curr-l))]
         [curr-z (/ (- curr-α (* prev-h prev-z)) curr-l)]
         [_ (displayln prev-h)]
         [next-c-l-z 
           (tridiagonal fn (cdr Xs) 
                        next-x curr-x 
                        curr-h
                        curr-μ curr-z
                        next-a curr-a
                        coefficients-table)]
         [_ (displayln next-c-l-z)]
         [next-c (first next-c-l-z)]
         [curr-c (- curr-z (* curr-μ next-c))]
         [curr-b (- (/ (- next-a curr-a) curr-h)
                    (/ (* curr-h (+ next-c (* 2 curr-c))) 3))]
         [curr-d (/ (- next-c curr-c) (* 3 curr-h))]
         [_ (hash-set! coefficients-table curr-x 
                       (list curr-a curr-b curr-c curr-d))])
        (list curr-c curr-l curr-z)))))
         
(define (alpha f x0 x1 x2)
  ( - (/ (* 3 (- (f x2) (f x1))) (- x2 x1))
      (/ (* 3 (- (f x1) (f x0))) (- x1 x0))))
