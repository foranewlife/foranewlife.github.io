# 个人学术网站

这是我的个人学术网站，专注于**计算生物学**和**密度泛函理论**的研究笔记和学术分享。

## 网站内容

- **学术笔记**：DFT理论、计算生物学方法、编程技巧
- **研究方向**：个人研究兴趣和项目介绍
- **学术简历**：教育背景、发表论文、会议报告
- **个人介绍**：联系方式和基本信息

## 技术栈

- **静态网站生成器**：Jekyll
- **主题**：Minima（简洁学术风格）
- **数学公式**：MathJax支持
- **部署平台**：GitHub Pages
- **自动化**：GitHub Actions

## 功能特色

✅ **数学公式支持**：完整的LaTeX数学公式渲染  
✅ **代码高亮**：多种编程语言语法高亮  
✅ **响应式设计**：适配桌面和移动设备  
✅ **SEO优化**：搜索引擎友好  
✅ **快速加载**：静态网站，加载速度快

## 本地开发

```bash
# 克隆仓库
git clone https://github.com/foranewlife/foranewlife.github.io.git
cd foranewlife.github.io

# 安装依赖
bundle install

# 本地运行
bundle exec jekyll serve

# 浏览器访问 http://localhost:4000
```

## 内容结构

```
.
├── _posts/           # 博客文章
├── _config.yml       # 网站配置
├── index.md          # 首页
├── about.md          # 关于页面
├── research.md       # 研究方向
├── notes.md          # 笔记索引
├── cv.md             # 学术简历
└── assets/           # 静态资源
```

## 写作指南

### 文章格式

每篇文章需要包含YAML Front Matter：

```yaml
---
layout: post
title: "文章标题"
date: YYYY-MM-DD
categories: [分类1, 分类2]
tags: [标签1, 标签2, 标签3]
author: "作者名"
---
```

### 数学公式

- 行内公式：`$E = mc^2$`
- 块级公式：
```latex
$$
\hat{H}\psi = E\psi
$$
```

### 代码块

```python
# Python代码示例
import numpy as np
print("Hello, World!")
```

## 部署说明

本网站使用GitHub Actions自动部署到GitHub Pages：

1. 推送代码到master分支
2. GitHub Actions自动构建Jekyll网站
3. 部署到 https://foranewlife.github.io

## 自定义配置

主要配置文件：

- `_config.yml`：网站基本设置
- `Gemfile`：Ruby依赖管理
- `.github/workflows/jekyll.yml`：自动部署配置

## 联系方式

如有问题或建议，欢迎通过以下方式联系：

- **Email**：your-email@example.com
- **GitHub Issues**：[提交问题](https://github.com/foranewlife/foranewlife.github.io/issues)

## 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件。

---

*最后更新：2024年1月*