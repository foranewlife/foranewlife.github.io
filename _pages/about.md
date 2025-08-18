---
permalink: /
title: "泳烁的个人网站"
author_profile: true
redirect_from: 
  - /about/
  - /about.html
---

欢迎来到我的个人网站！这里记录了我在研究中的学习笔记、技术分享和学术心得。

如有学术交流或合作意向，欢迎通过侧边栏的联系方式与我交流。

---

## 文章总览
### · [分类](/categories/)

{% comment %}
创建一个包含数量和名称的数组，用于排序
格式：[数量(补零到3位), 名称]，这样可以按字符串排序来实现按数量排序
{% endcomment %}
{% assign category_list = '' | split: '' %}
{% for category in site.categories %}
  {% assign count_padded = category[1].size | prepend: '000' | slice: -3, 3 %}
  {% assign item = count_padded | append: '|' | append: category[0] %}
  {% assign category_list = category_list | push: item %}
{% endfor %}
{% assign category_list = category_list | sort | reverse %}

<div class="category-stats">
  {% for item in category_list %}
    {% assign parts = item | split: '|' %}
    {% assign count = parts[0] | plus: 0 %}
    {% assign name = parts[1] %}
    <span class="category-item">
      <a href="/categories/#{{ name | slugify }}">{{ name }}</a>
      <span class="count">({{ count }})</span>
    </span>
    {% unless forloop.last %} • {% endunless %}
  {% endfor %}
</div>

### · [标签](/tags/)

{% comment %}
同样的方法处理标签
{% endcomment %}
{% assign tag_list = '' | split: '' %}
{% for tag in site.tags %}
  {% assign count_padded = tag[1].size | prepend: '000' | slice: -3, 3 %}
  {% assign item = count_padded | append: '|' | append: tag[0] %}
  {% assign tag_list = tag_list | push: item %}
{% endfor %}
{% assign tag_list = tag_list | sort | reverse %}

<div class="tag-stats">
  {% for item in tag_list %}
    {% assign parts = item | split: '|' %}
    {% assign count = parts[0] | plus: 0 %}
    {% assign name = parts[1] %}
    <span class="tag-item">
      <a href="/tags/#{{ name | slugify }}">{{ name }}</a>
      <span class="count">({{ count }})</span>
    </span>
    {% unless forloop.last %} • {% endunless %}
  {% endfor %}
</div>
---
{% if site.posts.size > 0 %}
<div class="page__latest-posts">
  <h2 class="page__latest-title">最新文章</h2>
  <div class="grid__wrapper">
    {% for post in site.posts limit:4 %}
      {% unless post.tags contains "草稿" %}
        {% include archive-single.html type="grid" %}
      {% endunless %}
    {% endfor %}
  </div>
</div>
{% endif %}