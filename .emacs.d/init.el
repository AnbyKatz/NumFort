(package-initialize)
(server-start)

(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/"))
(add-to-list 'package-archives '("gnu" . "https://elpa.gnu.org/packages/"))


;; Other files
(add-to-list 'custom-theme-load-path "~/.emacs.d/themes/")
(setq abbrev-file-name
      "~/.abbrev_defs")
(if (file-exists-p abbrev-file-name)
    (quietly-read-abbrev-file))


;; GUI-specific settings
(when (display-graphic-p)
  (load-theme 'abell-dark t)
  (add-to-list 'default-frame-alist '(fullscreen . maximized))
  (global-linum-mode t)
  (helm-mode 1)
  )

(set-language-environment "UTF-8")
(set-default-coding-systems 'utf-8)

;; testing

;; (defun ff ()
;;   "Prompt user to enter a string, with input history support."
;;   (interactive)
;;   (message "String is %s" (read-string "Enter your name:")))
;;  (insert "\\(\\:  \\:\\)")

(defun ask-name (x)
  "Ask name."
  (interactive "sEnter your name: ")
  (message "Name: %s" x))


;; ------------------------------------------------------------------- ;;

;; ---------f90 stuff------------------
(autoload 'f90-mode "f90" "Fortran 90 mode" t)

(defun my/f90-newline-yank ()
  (interactive)
  (newline)
  (yank)
  )

(defun my/f90-comment-header (title)
  "Inserts a commented title for f90"
  (interactive "sEnter a title: ")
  (defvar dash-len 1)
  (setq dash-len (/ (- 68 (length title)) 2))
  (indent-for-tab-command)
  (insert "! ")
  (dotimes (ii dash-len)
    (insert "-"))
  (if (= (mod (length title) 2) 1)
      (insert "-")
    )
  (insert title)
  (dotimes (ii dash-len)
    (insert "-"))
  )

(defun my/f90-comment-header-block (title)
  "Inserts a commented title block for f90"
  (interactive "sEnter a title: ")
  (defvar blank-len 1)
  (setq blank-len (/ (- 69 (length title)) 2))
  (newline)
  (dotimes (jj 5)
    (case jj
      ((0 4)
       (indent-for-tab-command)
       (insert "!")
       (dotimes (ii 69) (insert "-"))
       (insert "!")
       (newline))
      ((1 3)
       (indent-for-tab-command)
       (insert "!")
       (dotimes (ii 69) (insert " "))
       (insert "!")
       (newline))
      (2
       (indent-for-tab-command)
       (insert "!")
       (dotimes (ii blank-len)
	 (insert " "))
       (if (= (mod (length title) 2) 0)
	   (insert " ")
	 )
       (insert title)
       (dotimes (ii blank-len)
	 (insert " "))
       (insert "!")
       (newline))
  )))



;; (autoload 'f90-mode "f08" "Fortran 90 mode" t)
(defun my-f90-mode-hook ()
  (setq f90-font-lock-keywords f90-font-lock-keywords-3)
  '(f90-comment-region "!!!$")
  '(f90-indented-comment-re "!")
  (abbrev-mode 1)                       ; turn on abbreviation mode
  (turn-on-font-lock)                   ; syntax highlighting
  (auto-fill-mode 0)                    ; turn off auto-filling
  (local-set-key (kbd "<C-return>") 'my/f90-newline-yank)
  (local-set-key (kbd "H-t") 'my/f90-comment-header)
  (local-set-key (kbd "C-H-t") 'my/f90-comment-header-block)
  )
(add-hook 'f90-mode-hook 'my-f90-mode-hook)

;; ------------Fixed format Fortran stuff-----------------
;; remove the ^M at the end of lines
;; (defun remove-dos-eol ()
;;   "Do not show ^M in files containing mixed UNIX and DOS line endings."
;;   (interactive)
;;   (setq buffer-display-table (make-display-table))
;;   (aset buffer-display-table ?\^M []))
;; (add-hook 'fortran-mode-hook 'remove-dos-eol)

;; (defun my-fortran-mode-hook ()
;;   '(remove-dos-eol)
;;   )
;; (add-hook 'fortran-mode-hook 'my-fortran-mode-hook)


;; ------------Org-mode stuff------------
(autoload 'org-mode "org" "Org Mode" t)
(defun my-org-mode-hook ()
  (setq org-log-done t)
  (define-key global-map "\C-cl" 'org-store-link)
  (define-key global-map "\C-ca" 'org-agenda)
  (visual-line-mode 1)
  (org-indent-mode 1)
  (abbrev-mode 1)
  (setq org-src-fontify-natively t)
  ;; (set-language-environment "UTF-8")
  ;; (set-default-coding-systems 'utf-8)
  )
(add-hook 'org-mode-hook 'my-org-mode-hook)

(setq org-src-fontify-natively t)

;; ------------pdf tools------------
(add-to-list 'auto-mode-alist '("\\.pdf\\'" . pdf-view-mode))
;; make midnight mode colours nice
(setq pdf-view-midnight-colors (cons (face-foreground 'default) (face-background 'default)))
(defun my-pdf-view-mode-hook ()
  (pdf-view-midnight-minor-mode 1)
  (linum-mode 0)
  )
(add-hook 'pdf-view-mode-hook 'my-pdf-view-mode-hook)


(pdf-tools-install)
;; to use pdfview with auctex
(setq TeX-view-program-selection '((output-pdf "PDF Tools"))
      TeX-view-program-list '(("PDF Tools" TeX-pdf-tools-sync-view))
      TeX-source-correlate-start-server t) ;; not sure if last line is neccessary

;; to have the buffer refresh after compilation
(add-hook 'TeX-after-compilation-finished-functions
	  #'TeX-revert-document-buffer)



;; --------------------TeX--------------------
(defun TeX-inline-math-abell()
  (interactive)
  (insert "\\(\\:  \\:\\)")
  (backward-char 5))
(defun TeX-fullline-math-abell()
  (interactive)
  (insert "\\[  \\]")
  (backward-char 3))

(defun TeX-align-newline-abell()
  (interactive)
  (insert "\\\\")
  (newline)
  (insert "&= ")
  (indent-for-tab-command))

(require 'tex)
(defun my-LaTeX-mode-hook ()
  (setq TeX-auto-save t)
  (setq TeX-parse-self t)
  (setq-default TeX-master nil)
  (setq TeX-PDF-mode t)
  (visual-line-mode 1)
  (flyspell-mode 1)
  (LaTeX-math-mode 1)
  (TeX-source-correlate-mode 1)
  (outline-minor-mode 1)
  (rainbow-delimiters-mode 1)
  (local-set-key (kbd "C-c m") 'TeX-inline-math-abell)
  (local-set-key (kbd "C-c m") 'TeX-inline-math-abell)
  (local-set-key (kbd "C-c M-m") 'TeX-fullline-math-abell)
  (local-set-key (kbd "C-M-=") '(lambda () (interactive) (insert "&= ")))
  (local-set-key (kbd "C-c b") 'tex-latex-block)
  (local-set-key (kbd "<C-tab>") 'outline-toggle-children)
  (local-set-key (kbd "<C-return>") 'TeX-align-newline-abell)
  )
(add-hook 'LaTeX-mode-hook 'my-LaTeX-mode-hook)


;; ------------Emacs-Lisp------------
(define-key emacs-lisp-mode-map (kbd "C-c C-a") 'eval-buffer)
(define-key emacs-lisp-mode-map (kbd "C-c C-r") 'eval-region)

;; ------------Wolfram Alpha---------------
(setq wolfram-alpha-app-id 'H8TGH2-5LXPKP5V9J)

;; ------------Helm Spotify Plus-------------
(require 'helm-spotify-plus)
(global-set-key (kbd "H-s H-SPC") 'helm-spotify-plus-toggle-play-pause)
(global-set-key (kbd "H-s H-n") 'helm-spotify-plus-next)
(global-set-key (kbd "H-s H-p") 'helm-spotify-plus-previous)

;; ---------------Simpleclip---------------
(simpleclip-mode 1)
(require 'simpleclip)
(global-set-key (kbd "H-c") 'simpleclip-copy)
(global-set-key (kbd "H-x") 'simpleclip-cut)
(global-set-key (kbd "H-v") 'simpleclip-paste)

;; ------------------------------------------------------------------- ;;

;; Use sensible-defaults.el for some basic settings.
(load-file "~/.emacs.d/code/sensible-defaults.el")
(sensible-defaults/use-all-settings)
(sensible-defaults/use-all-keybindings)
(sensible-defaults/backup-to-temp-directory)

;; fuck tabs
(setq-default indent-tabs-mode 0)

;; Convert a file with DOS line endings to a UNIX file
(defun dos2unix (buffer)
  "Automate M-% C-q C-m RET C-q C-j RET"
  (interactive "*b")
  (save-excursion
    (goto-char (point-min))
    (while (search-forward (string ?\C-m) nil t)
      (replace-match (string ?\C-j) nil t))))


;; ;; folding minor mode
;; (load-file "~/.emacs.d/code/folding.el")
;; (load "folding" 'nomessage 'noerror)
;; (folding-mode-add-find-file-hook)
;; (folding-add-to-marks-list 'ruby-mode "#{{{" "#}}}" nil t)
;; (folding-add-to-marks-list 'f90-mode "!{" "!}" nil t)


;; ----------------Ace Window------------------
(global-set-key (kbd "H-o") 'ace-select-window)
;; macros for quick switch - disgusting
(fset 'ace-switch-1
   [?\H-o ?1])
(fset 'ace-switch-2
   [?\H-o ?2])
(fset 'ace-switch-3
   [?\H-o ?3])
(fset 'ace-switch-4
   [?\H-o ?4])
(global-set-key (kbd "H-1") 'ace-switch-1)
(global-set-key (kbd "H-2") 'ace-switch-2)
(global-set-key (kbd "H-3") 'ace-switch-3)
(global-set-key (kbd "H-4") 'ace-switch-4)

;; (defun my/ace-switch-to-window (win_num)
;;   (interactive)
;;   (aw-switch-to-window win_num)
;;   )
;; aw-switch-to-window(1)


;(transient-mark-mode 0)

;; Delete selection upon backspacing or typing.
(delete-selection-mode 0)

;; Don't show menu or scroll bar
(tool-bar-mode 0)
(menu-bar-mode 0)
(when window-system
  (scroll-bar-mode -1))


;; (setq split-height-threshold nil)
;; (setq split-width-threshold 0)

;; Open GUI in full screen by default
;; (set-frame-parameter nil 'fullscreen 'fullboth)

;; Set vertical bar cursor
;(setq-default cursor-type 'bar)


;; Emacs has a silly warning bell by default. This gets rid of it.
(setq ring-bell-function 'ignore)


;; Fix Emacs' mouse scrolling behaviour
(setq scroll-conservatively 100) ;; When cursor moves outside window, don't jump erratically
(setq mouse-wheel-scroll-amount '(1 ((shift) . 1))) ;; one line at a time
;; (setq mouse-wheel-progressive-speed nil) ;; don't accelerate scrolling
(setq mouse-wheel-follow-mouse 't) ;; scroll window under mouse

;; Font options
;(set-frame-font "Menlo 14")

;; Some handy functions to have
(defun my/view-buffer-name ()
  "Display the filename of the current buffer."
  (interactive)
  (message (buffer-file-name)))
(global-set-key (kbd "H-b") 'my/view-buffer-name)


(defun my/rename-file (new-name)
  (interactive "FNew name: ")
  (let ((filename (buffer-file-name)))
    (if filename
        (progn
          (when (buffer-modified-p)
            (save-buffer))
          (rename-file filename new-name t)
          (kill-buffer (current-buffer))
          (find-file new-name)
          (message "Renamed '%s' -> '%s'" filename new-name))
      (message "Buffer '%s' isn't backed by a file!" (buffer-name)))))


;; Highlight the current line in GUI
;(when window-system
;  (global-hl-line-mode))


;; Display time in mode line
(setq display-time-string-forms
      '((propertize (format-time-string " %b %d, %l:%M%P" now) 'face 'bold)))
(setq display-time-and-date t)
(display-time-mode 1)

;; Abbreviate all 'Yes/No' prompts to 'y/n'
(fset 'yes-or-no-p 'y-or-n-p)

;; Linum is great, but it slows certain modes down. linum-off redefines
;; global-linum-mode to exclude modes
;(require 'linum-off)
(require 'linum)
(setq linum-disabled-modes-list '(eshell-mode wl-summary-mode
     compilation-mode text-mode dired-mode pdf-view-mode
     doc-view-mode shell-mode pdf-view-mode image-mode))

;; Fix startup behaviour. Don't show startup screen, replace with dashboard.
(setq inhibit-startup-screen t)

;; Enable quick access to emacs init file with "C-c e"
(defun my/visit-emacs-config ()
  (interactive)
  (find-file "~/.emacs.d/init.el"))
;; Open emacs init file in other window with C-c M-e
(defun my/visit-emacs-config-other-window ()
  (interactive)
  (find-file-other-window "~/.emacs.d/init.el"))
(global-set-key (kbd "C-c e") 'my/visit-emacs-config)
(global-set-key (kbd "C-c M-e") 'my/visit-emacs-config-other-window)

;; Settings changed through Emacs interface are stored in separate file
(setq custom-file "~/.emacs.d/custom.el")
(load custom-file :noerror)

;; camelCase recognition
(global-subword-mode)

;; c global tab width
(setq tab-width 5)

;; Avoid truncation of emacs term (default 2048)
;; (setq 'term-buffer-maximum-size 0)

;; Move lines function
(defun move-line-up ()
  (interactive)
  (transpose-lines 1)
  (previous-line 2))

(defun move-line-down ()
  (interactive)
  (forward-line 1)
  (transpose-lines 1)
  (previous-line 1))

(global-set-key (kbd "M-<up>") 'move-line-up)
(global-set-key (kbd "M-<down>") 'move-line-down)

;; Helm
(require 'helm)
(require 'helm-config)
(global-set-key (kbd "M-x") 'helm-M-x)
(global-set-key (kbd "C-x C-f") 'helm-find-files)
(global-set-key (kbd "C-x C-b") 'helm-buffers-list)

;; Enable disabled commands
(put 'narrow-to-region 'disabled nil)

;; duplicate line
(defun duplicate-current-line-or-region (arg)
  "Duplicates the current line or region ARG times.
If there's no region, the current line will be duplicated. However, if
there's a region, all lines that region covers will be duplicated."
  (interactive "p")
  (let (beg end (origin (point)))
    (if (and mark-active (> (point) (mark)))
        (exchange-point-and-mark))
    (setq beg (line-beginning-position))
    (if mark-active
        (exchange-point-and-mark))
    (setq end (line-end-position))
    (let ((region (buffer-substring-no-properties beg end)))
      (dotimes (i arg)
        (goto-char end)
        (newline)
        (insert region)
        (setq end (point)))
      (goto-char (+ origin (* (length region) arg) arg)))))
(global-set-key (kbd "H-d") 'duplicate-current-line-or-region)


(defun kill-buffer-and-frame ()
  (interactive)
  (kill-this-buffer)
  (delete-frame))


(defun my/repeat-last-shell-command()
  (interactive)
  (shell-command (cadr (assoc 'shell-command command-history))))


(defun my/cat-this-file()
  (interactive)
  (defvar thisfile buffer-file-name)
  (shell-command (concat "cat " thisfile)))


;; custom keybindings
(global-set-key (kbd "M-n") 'forward-paragraph)
(global-set-key (kbd "M-p") 'backward-paragraph)
(global-set-key (kbd "M-]") 'other-frame)
(global-set-key (kbd "M-[") 'other-window)
(global-set-key (kbd "C-x 4 k") 'kill-buffer-and-window)
(global-set-key (kbd "C-x 5 k") 'kill-buffer-and-frame)
(global-set-key (kbd "H-M-z") (lambda () (interactive) (shell-command "./macro_z")))
(global-set-key (kbd "H-M-c") (lambda () (interactive) (shell-command "./macro_c")))
(global-set-key (kbd "H-M-i") (lambda () (interactive) (shell-command "./macro_i")))
(global-set-key (kbd "H-M-f") (lambda () (interactive) (shell-command "./macro_f")))
(global-set-key (kbd "H-M-x") (lambda () (interactive) (shell-command "./macro_x")))
(global-set-key (kbd "H-M-p") (lambda () (interactive) (shell-command "./macro_p")))
(global-set-key (kbd "<menu>") 'shell-command)
(global-set-key (kbd "C-H-<menu>") 'my/repeat-last-shell-command)
