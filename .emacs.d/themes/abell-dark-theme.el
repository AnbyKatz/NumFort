(deftheme abell-dark
  "abell-dark theme")

(custom-theme-set-faces
 'abell-dark

 ;;; Normal Text
 '(default ((t (:background "black" :foreground "White"))))
 '(mouse ((t (:foreground "grey85"))))
 '(cursor ((t (:background "grey85"))))

 ;;; Comment colour
 '(font-lock-comment-face ((t (:italic t :foreground "grey60"))))
 ;;; String colour
 '(font-lock-string-face ((t (:foreground "Magenta"))))
 ;;; call, write, allocate, open, do, if etc
 '(font-lock-keyword-face ((t (:bold t :foreground "Cyan"))))
 ;;; '(font-lock-keyword-face
 '(font-lock-warning-face ((t (:bold t :foreground "Pink"))))
 ;;; No use in fortran?
 '(font-lock-constant-face ((t (:foreground "DarkGreen"))))
 ;;; Variable types (real, int etc) ForestGreen
 '(font-lock-type-face ((t (:foreground "green3"))))
 ;;; Variable names
 '(font-lock-variable-name-face ((t (:foreground "goldenrod1"))));;;"DarkGoldenrod"))))
 ;;; Program name, modules, functions etc
 '(font-lock-function-name-face ((t (:foreground "SpringGreen"))))
 ;; '(font-lock-function-name-face ((t (:foreground "blue2"))))
 ;;; other stuff i couldnt be bothered figuring out lol
 '(font-lock-builtin-face ((t (:foreground "SkyBlue"))))
 '(highline-face ((t (:background "grey12"))))
 '(setnu-line-number-face ((t (:background "Grey15" :foreground "White" :bold t))))
 '(show-paren-match-face ((t (:background "grey30"))))
 '(region ((t (:background "grey15"))))
 '(highlight ((t (:background "blue"))))
 '(secondary-selection ((t (:background "navy"))))
 '(widget-field-face ((t (:background "navy"))))
 '(widget-single-line-field-face ((t (:background "royalblue")))))

;;;###autoload
(when load-file-name
  (add-to-list 'custom-theme-load-path
               (file-name-as-directory (file-name-directory load-file-name))))

(provide-theme 'abell-dark)
