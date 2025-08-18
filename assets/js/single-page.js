/**
 * Single page functionality
 * Includes back to top button and table of contents functionality
 */

document.addEventListener('DOMContentLoaded', function() {
  // Back to top button functionality
  const backToTopBtn = document.getElementById('back-to-top');
  
  if (backToTopBtn) {
    // Show/hide button based on scroll position
    function toggleBackToTopBtn() {
      if (window.pageYOffset > 300) {
        backToTopBtn.classList.add('show');
      } else {
        backToTopBtn.classList.remove('show');
      }
    }
    
    // Smooth scroll to top
    function scrollToTop() {
      window.scrollTo({
        top: 0,
        behavior: 'smooth'
      });
    }
    
    // Listen for scroll events
    window.addEventListener('scroll', toggleBackToTopBtn);
    
    // Listen for click events
    backToTopBtn.addEventListener('click', scrollToTop);
  }
  
  // Table of contents functionality
  if (document.querySelector('.toc-sidebar')) {
    // Generate table of contents
    function generateTOC() {
      const content = document.querySelector('.page__content');
      const tocContent = document.getElementById('toc-content');
      
      if (!content || !tocContent) return;
      
      const headings = content.querySelectorAll('h1, h2, h3, h4, h5, h6');
      if (headings.length === 0) {
        document.querySelector('.toc-sidebar').style.display = 'none';
        return;
      }
      
      // Build hierarchical structure
      const tocTree = buildTOCTree(headings);
      tocTreeGlobal = tocTree; // Save global reference
      const tocHTML = generateTOCHTML(tocTree);
      
      tocContent.innerHTML = tocHTML;
    }
    
    // Build TOC tree structure
    function buildTOCTree(headings) {
      const tree = [];
      const stack = [];
      
      headings.forEach((heading, index) => {
        // Add ID to headings if they don't have one
        if (!heading.id) {
          const text = heading.textContent.trim();
          heading.id = `heading-${index}-${text.replace(/[^\u4e00-\u9fa5a-zA-Z0-9]/g, '').substring(0, 20)}`;
        }
        
        const level = parseInt(heading.tagName.substring(1));
        const text = heading.textContent.trim();
        
        const node = {
          id: heading.id,
          text: text,
          level: level,
          children: []
        };
        
        // Find appropriate parent node
        while (stack.length > 0 && stack[stack.length - 1].level >= level) {
          stack.pop();
        }
        
        if (stack.length === 0) {
          tree.push(node);
        } else {
          stack[stack.length - 1].children.push(node);
        }
        
        stack.push(node);
      });
      
      return tree;
    }
    
    // Global variable to store TOC tree for parent highlighting
    let tocTreeGlobal = null;
    
    // Find all parent IDs of a node
    function findParentIds(nodeId, tree, parents = []) {
      for (let node of tree) {
        if (node.id === nodeId) {
          return parents;
        }
        if (node.children && node.children.length > 0) {
          const result = findParentIds(nodeId, node.children, [...parents, node.id]);
          if (result !== null) {
            return result;
          }
        }
      }
      return null;
    }
    
    // Highlight parent links
    function highlightParents(currentId) {
      if (!tocTreeGlobal) return;
      
      const parentIds = findParentIds(currentId, tocTreeGlobal);
      if (parentIds) {
        parentIds.forEach(parentId => {
          const parentLink = document.querySelector(`.toc-link[href="#${parentId}"]`);
          if (parentLink) {
            parentLink.classList.add('parent-active');
          }
        });
      }
    }
    
    // Scroll active TOC item into view
    function scrollTocToActiveItem(activeLink) {
      const tocNav = document.querySelector('.toc-nav');
      if (!tocNav || !activeLink) return;
      
      const tocNavRect = tocNav.getBoundingClientRect();
      const activeLinkRect = activeLink.getBoundingClientRect();
      
      // Calculate position relative to toc-nav container
      const relativeTop = activeLinkRect.top - tocNavRect.top + tocNav.scrollTop;
      const relativeBottom = activeLinkRect.bottom - tocNavRect.top + tocNav.scrollTop;
      
      // Check if scrolling is needed
      const containerHeight = tocNav.clientHeight;
      const scrollTop = tocNav.scrollTop;
      const visibleTop = scrollTop;
      const visibleBottom = scrollTop + containerHeight;
      
      // If active item is above visible area
      if (relativeTop < visibleTop) {
        tocNav.scrollTo({
          top: relativeTop - 20, // Add some margin
          behavior: 'smooth'
        });
      }
      // If active item is below visible area
      else if (relativeBottom > visibleBottom) {
        tocNav.scrollTo({
          top: relativeBottom - containerHeight + 20, // Add some margin
          behavior: 'smooth'
        });
      }
    }
    
    // Generate TOC HTML
    function generateTOCHTML(tree) {
      if (!tree || tree.length === 0) return '';
      
      let html = '<ul class="toc-list">';
      
      tree.forEach(node => {
        const hasChildren = node.children && node.children.length > 0;
        
        html += `<li class="toc-level-${node.level} ${hasChildren ? 'toc-has-children' : ''}">`;
        
        html += `<a href="#${node.id}" class="toc-link">`;
        
        if (hasChildren) {
          html += `<span class="toc-toggle" data-target="toc-${node.id}">
            <i class="fa fa-chevron-down"></i>
          </span>`;
        }
        
        html += `${node.text}</a>`;
        
        if (hasChildren) {
          html += `<div class="toc-children" id="toc-${node.id}">
            ${generateTOCHTML(node.children)}
          </div>`;
        }
        
        html += '</li>';
      });
      
      html += '</ul>';
      return html;
    }
    
    // Update TOC highlighting based on scroll
    function updateTOCHighlight() {
      const headings = document.querySelectorAll('.page__content h1[id], .page__content h2[id], .page__content h3[id], .page__content h4[id], .page__content h5[id], .page__content h6[id]');
      const tocLinks = document.querySelectorAll('.toc-link');
      
      let currentHeading = null;
      const scrollTop = window.pageYOffset + 100;
      
      // Find current visible heading
      for (let i = headings.length - 1; i >= 0; i--) {
        if (headings[i].offsetTop <= scrollTop) {
          currentHeading = headings[i];
          break;
        }
      }
      
      // Update highlighting
      tocLinks.forEach(link => {
        link.classList.remove('active');
        link.classList.remove('parent-active');
      });
      
      if (currentHeading) {
        const activeLink = document.querySelector(`.toc-link[href="#${currentHeading.id}"]`);
        if (activeLink) {
          activeLink.classList.add('active');
          
          // Highlight parent links
          highlightParents(currentHeading.id);
          
          // Scroll TOC to visible area
          scrollTocToActiveItem(activeLink);
        }
      }
    }
    
    // Handle toggle collapse/expand
    function handleTOCToggle(e) {
      if (e.target.closest('.toc-toggle')) {
        e.preventDefault();
        const toggle = e.target.closest('.toc-toggle');
        const targetId = toggle.getAttribute('data-target');
        const children = document.getElementById(targetId);
        const icon = toggle.querySelector('i');
        
        if (children) {
          children.classList.toggle('collapsed');
          icon.classList.toggle('fa-chevron-down');
          icon.classList.toggle('fa-chevron-right');
        }
      }
    }
    
    // Handle smooth scroll to headings
    function handleTOCClick(e) {
      if (e.target.classList.contains('toc-link')) {
        e.preventDefault();
        const targetId = e.target.getAttribute('href').substring(1);
        const targetElement = document.getElementById(targetId);
        if (targetElement) {
          window.scrollTo({
            top: targetElement.offsetTop - 80,
            behavior: 'smooth'
          });
          
          // Update highlighting and TOC scroll immediately
          setTimeout(() => {
            updateTOCHighlight();
          }, 100);
        }
      }
    }
    
    // Mobile TOC toggle functionality
    function initMobileTOC() {
      const mobileToggle = document.getElementById('mobile-toc-toggle');
      const tocSidebar = document.querySelector('.toc-sidebar');
      
      // Auto-configure top position
      function updateTocPosition() {
        const masthead = document.querySelector('.masthead');
        if (masthead && tocSidebar) {
          const mastheadHeight = masthead.offsetHeight;
          tocSidebar.style.top = mastheadHeight + 'px';
          tocSidebar.style.height = `calc(100vh - ${mastheadHeight}px)`;
        }
      }
      
      // Update position after page load and on window resize
      setTimeout(updateTocPosition, 100); // Ensure DOM is fully loaded
      window.addEventListener('resize', updateTocPosition);
      
      if (mobileToggle && tocSidebar) {
        // Toggle TOC display on button click
        mobileToggle.addEventListener('click', function() {
          tocSidebar.classList.toggle('show');
          
          // Change icon
          const icon = mobileToggle.querySelector('i');
          if (tocSidebar.classList.contains('show')) {
            icon.className = 'fa fa-times';
          } else {
            icon.className = 'fa fa-list';
          }
        });
        
        // Auto-close mobile TOC after clicking a link
        tocSidebar.addEventListener('click', function(e) {
          if (e.target.classList.contains('toc-link') && window.innerWidth <= 768) {
            tocSidebar.classList.remove('show');
            const icon = mobileToggle.querySelector('i');
            icon.className = 'fa fa-list';
          }
        });
        
        // Close TOC when clicking outside
        document.addEventListener('click', function(e) {
          if (window.innerWidth <= 768 && 
              tocSidebar.classList.contains('show') && 
              !tocSidebar.contains(e.target) && 
              !mobileToggle.contains(e.target)) {
            tocSidebar.classList.remove('show');
            const icon = mobileToggle.querySelector('i');
            icon.className = 'fa fa-list';
          }
        });
      }
    }
    
    // Initialize TOC functionality
    generateTOC();
    initMobileTOC();
    window.addEventListener('scroll', updateTOCHighlight, { passive: true });
    
    // Bind click events (collapse and navigation)
    const tocSidebar = document.querySelector('.toc-sidebar');
    if (tocSidebar) {
      tocSidebar.addEventListener('click', handleTOCToggle);
      tocSidebar.addEventListener('click', handleTOCClick);
    }
    
    // Initial highlighting
    setTimeout(updateTOCHighlight, 100);
  }
});