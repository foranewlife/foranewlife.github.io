/**
 * Archive Scroll Spy - 归档页面滚动监听脚本
 * 监听页面滚动，高亮当前可见章节对应的侧边栏链接
 */

(function() {
  'use strict';

  // 配置选项
  const options = {
    offset: 100, // 距离顶部的偏移量
    throttleDelay: 100 // 节流延迟时间
  };

  let isThrottled = false;

  // 节流函数
  function throttle(func, delay) {
    return function(...args) {
      if (!isThrottled) {
        func.apply(this, args);
        isThrottled = true;
        setTimeout(() => {
          isThrottled = false;
        }, delay);
      }
    };
  }

  // 获取当前可见的标题元素
  function getCurrentVisibleSection() {
    const sections = document.querySelectorAll('.archive__subtitle');
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    
    let currentSection = null;
    let lastPassedSection = null;

    // 从上到下遍历所有标题
    for (let i = 0; i < sections.length; i++) {
      const section = sections[i];
      const rect = section.getBoundingClientRect();
      const sectionTop = rect.top + scrollTop;
      
      // 如果标题还在视窗上方（还没滚动到）
      if (sectionTop > scrollTop + options.offset) {
        // 如果是第一个标题且在视窗上方，选择它
        if (i === 0) {
          currentSection = section;
        }
        // 否则选择上一个已经过了的标题
        else if (lastPassedSection) {
          currentSection = lastPassedSection;
        }
        break;
      }
      // 标题已经被滚动过了，记录它
      else {
        lastPassedSection = section;
        currentSection = section;
      }
    }

    // 如果所有标题都被滚动过了，选择最后一个
    if (!currentSection && sections.length > 0) {
      currentSection = sections[sections.length - 1];
    }

    return currentSection;
  }

  // 更新侧边栏高亮状态
  function updateActiveLink(activeSection) {
    // 移除所有active类
    const allLinks = document.querySelectorAll('.category-nav a');
    allLinks.forEach(link => link.classList.remove('active'));

    if (activeSection) {
      const sectionId = activeSection.id;
      const activeLink = document.querySelector(`.category-nav a[href="#${sectionId}"]`);
      
      if (activeLink) {
        activeLink.classList.add('active');
        
        // 可选：滚动侧边栏确保活跃链接可见
        const sidebar = activeLink.closest('.archive-sidebar');
        if (sidebar) {
          const linkRect = activeLink.getBoundingClientRect();
          const sidebarRect = sidebar.getBoundingClientRect();
          
          if (linkRect.bottom > sidebarRect.bottom || linkRect.top < sidebarRect.top) {
            activeLink.scrollIntoView({
              behavior: 'smooth',
              block: 'center'
            });
          }
        }
      }
    }
  }

  // 处理滚动事件
  const handleScroll = throttle(() => {
    const currentSection = getCurrentVisibleSection();
    updateActiveLink(currentSection);
  }, options.throttleDelay);

  // 处理链接点击事件
  function handleLinkClick(event) {
    const link = event.target;
    if (link.matches('.category-nav a[href^="#"]')) {
      event.preventDefault();
      
      const targetId = link.getAttribute('href').substring(1);
      const targetElement = document.getElementById(targetId);
      
      if (targetElement) {
        const targetTop = targetElement.getBoundingClientRect().top + 
                         window.pageYOffset - options.offset;
        
        window.scrollTo({
          top: targetTop,
          behavior: 'smooth'
        });
        
        // 立即更新高亮状态
        setTimeout(() => {
          updateActiveLink(targetElement);
        }, 100);
      }
    }
  }

  // 初始化函数
  function init() {
    // 检查是否存在必要的元素
    const sidebar = document.querySelector('.archive-sidebar');
    const sections = document.querySelectorAll('.archive__subtitle');
    
    if (!sidebar || sections.length === 0) {
      return; // 如果不是归档页面，则不启用功能
    }

    // 绑定事件监听器
    window.addEventListener('scroll', handleScroll, { passive: true });
    window.addEventListener('resize', throttle(handleScroll, 250), { passive: true });
    sidebar.addEventListener('click', handleLinkClick);

    // 初始化高亮状态
    setTimeout(() => {
      const initialSection = getCurrentVisibleSection();
      updateActiveLink(initialSection);
    }, 100);

    console.log('Archive Scroll Spy initialized');
  }

  // DOM加载完成后初始化
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }

})();