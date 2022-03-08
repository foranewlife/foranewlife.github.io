# 陈泳烁的个人主页

本手册由Sphinx生成,并使用[Executable Book Project 主题](https://ebp.jupyterbook.org/),可以解析markdown与Latex公式。
[Sphinx官网](https://www.sphinx-doc.org/en/master/)、
[Sphinx中文手册](https://zh-sphinx-doc.readthedocs.io/en/latest/index.html)。

本手册参考项目为 [phonopy手册](http://phonopy.github.io/phonopy/)。

以下为与本项目相关的知识的最小合集

## Sphinx

### 安装Sphinx

```bash
pip install sphinx myst-parser sphinx-book-theme sphinxcontrib-mermaid
# or
pip install -r requirements.txt
#option
# pip install sphinx-math-dollar
# pip install -U sphinx-mathjax-offline
```

### 生成静态网页

```bash
make html
```

### Latex 相关

```bash
make latex
```

修改latex 相关参数

需要修改conf.py

参见

[www.sphinx-doc.org/en/master/latex.html](https://www.sphinx-doc.org/en/master/latex.html)
