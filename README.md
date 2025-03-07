# stark-anatomy

STARK tutorial with supporting code in python

Outline:

- introduction
- overview of STARKs
- basic tools -- algebra and polynomials
- FRI low degree test
- STARK information theoretical protocol
- speeding things up with NTT and preprocessing

Visit the Github Pages website here:
[https://aszepieniec.github.io/stark-anatomy/](https://aszepieniec.github.io/stark-anatomy/)

## Follow-up

Be sure to check out the
[next tutorial](https://github.com/aszepieniec/stark-brainfuck) where we
implement a STARK engine for a VM running Brainfuck. And our "real", functional,
practical ZK-STARK VM, [Triton VM](https://triton-vm.org/).

## Runner the tutorial locally

### Prerequisites

- Python 3.8 or higher
- [uv](https://github.com/astral-sh/uv)
- [trunk](https://docs.trunk.io/code-quality)

### Installation with uv (recommended)

1. Install uv and trunk if you don't have them already:

```bash
curl -sSf https://astral.sh/uv/install.sh | bash
curl https://get.trunk.io -fsSL | bash
```

2. Clone the repository:

```bash
git clone https://github.com/aszepieniec/stark-anatomy.git
cd stark-anatomy
```

3. Create and activate a virtual environment:

```bash
uv sync
```

### Running the tests

After installation, you can run the examples from the `src/stark_anatomy`
directory:

```bash
uv run pytest
```

## Running the website locally

1.  Install ruby
2.  Install bundler
3.  Change directory to `docs/` and install Jekyll: `$> sudo bundle install`
4.  Run Jekyll: `$> bundle exec jekyll serve`
5.  Surf to [http://127.0.0.1:4000/](http://127.0.0.1:4000/)

## LaTeX and Github Pages

Github-Pages uses Kramdown as the markdown processor. Kramdown does not support
LaTeX. Instead, there is a javascript header that loads MathJax, parses the
page, and replaces LaTeX maths instructions with their proper formulae. Here is
how to do it:

1. Open `_includes/head-custom.html` and paste the following code:

```javascript
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
    TeX: {
      equationNumbers: {
        autoNumber: "AMS"
      }
    },
    tex2jax: {
    inlineMath: [ ['$', '$'], ['\\(', '\\)'] ],
    displayMath: [ ['$$', '$$'], ['\\[', '\\]'] ],
    processEscapes: true,
  }
});
MathJax.Hub.Register.MessageHook("Math Processing Error",function (message) {
	  alert("Math Processing Error: "+message[1]);
	});
MathJax.Hub.Register.MessageHook("TeX Jax - parse error",function (message) {
	  alert("Math Processing Error: "+message[1]);
	});
</script>
<script type="text/javascript" async
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-MML-AM_CHTML">
</script>
```

Jekyll, the site engine used by Github Pages, will load this header
automatically. There is no need to change the `_config.yml` file.

Note that Kramdown interprets every underscore (`_`) that is followed by a
non-whitespace character, as starting an emphasized piece of text. This
interpretation interferences with subscript in LaTeX formulae, which also uses
underscores. The workaround is to re-write the LaTeX formulas by introducing a
space after every underscore. Also, consider replacing:

- `\{` by `\lbrace` and `\}` by `\rbrace`,
- `|` by `\vert`.
