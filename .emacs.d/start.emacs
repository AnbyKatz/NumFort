(package-initialize)

(require 'package)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/"))
(add-to-list 'package-archives '("gnu" . "https://elpa.gnu.org/packages/"))

;; Use sensible-defaults.el for some basic settings.
(load-file "~/.emacs.d/code/sensible-defaults.el")
(sensible-defaults/use-all-settings)
(sensible-defaults/use-all-keybindings)
(sensible-defaults/backup-to-temp-directory)

;; for extra emacs themes
(add-to-list 'custom-theme-load-path "~/.emacs.d/themes/")
;; abbrev-mode customisations
(setq abbrev-file-name
      "~/.abbrev_defs")
(if (file-exists-p abbrev-file-name)
    (quietly-read-abbrev-file))



;; --------------------f90 customisation--------------------
(defun my-f90-mode-hook ()
  (setq f90-font-lock-keywords f90-font-lock-keywords-3)
  '(f90-comment-region "!!")
  '(f90-indented-comment-re "!")
  (turn-on-font-lock)                   ; syntax highlighting
  (auto-fill-mode 0)                    ; turn off auto-filling
  )
(add-hook 'f90-mode-hook 'my-f90-mode-hook)


;; --------------------Misc Customisations--------------------

;; no startup screen
(inhibit-startup-screen t)

;; stop annoying cursor blinking
(blink-cursor-mode nil)

;; highlight matching brackets
(show-paren-mode t)

;; Don't show menu or scroll bar
(tool-bar-mode 0)
(menu-bar-mode 0)
(when window-system
  (scroll-bar-mode -1))

;; Fix Emacs' mouse scrolling behaviour
(setq scroll-conservatively 100) ;; When cursor moves outside window, don't jump erratically
(setq mouse-wheel-scroll-amount '(1 ((shift) . 1))) ;; one line at a time
(setq mouse-wheel-progressive-speed nil) ;; don't accelerate scrolling
(setq mouse-wheel-follow-mouse 't) ;; scroll window under mouse

;; Emacs has a silly warning bell by default. This gets rid of it.
(setq ring-bell-function 'ignore)

;; Settings changed through Emacs interface are stored in separate file
(setq custom-file "~/.emacs.d/custom.el")
(load custom-file :noerror)

;; camelCase recognition
(global-subword-mode)
