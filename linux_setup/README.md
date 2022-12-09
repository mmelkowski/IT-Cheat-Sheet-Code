# Linux Fedora OS: A quick setup guide

This folder describre my quick setup for fedora OS.

## Default fedora

Fedora [website](https://getfedora.org/fr/workstation/) to download the OS.

### visuel

* [dashtopanel](https://extensions.gnome.org/extension/1160/dash-to-panel/)
* gnome-tweaks (*sudo dnf install gnome-tweaks*)

## Custom Desktop env

See [spins](https://spins.fedoraproject.org/) for other environment for fedora such as:

* [KDE plasma](https://spins.fedoraproject.org/kde/)
* [LXQT](https://spins.fedoraproject.org/fr/lxqt/)

## Console: Oh My Zsh

### customization

Use of agnoster theme with added time inside the command line.

Change zsh theme in ~/.zshrc:

```bash
ZSH_THEME="agnoster_custom"
```

The theme can be found in the zsh_theme folder and be added to the `.oh-my-zsh/themes` folder.

It has an added time function:

```bash
# Context: Date
prompt_date() {
  prompt_segment black default "$(date +'%Y-%m-%d %H:%M:%S')"
}

## Main prompt
build_prompt() {
  RETVAL=$?
  prompt_status
  prompt_virtualenv
  prompt_aws
  prompt_context
  prompt_date
  prompt_dir
  prompt_git
  prompt_bzr
  prompt_hg
  prompt_end
}
```

## Standard alias

edit or create ~/.bash_aliases

```bash
alias ss='git status'
```

and source it in ~/.zshrc:

```bash
# source personnal alias
source ~/.bash_aliases
```