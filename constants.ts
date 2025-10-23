import { Project } from './types';

export const SOCIAL_LINKS = {
  // --- IMPORTANT: Replace these placeholder links with your own! ---
  linkedIn: 'https://www.linkedin.com/in/your-profile-here',
  googleScholar: 'https://scholar.google.com/citations?user=your-id-here',
  github: 'https://github.com/your-username-here',
  // Replace '#' with a link to your actual resume PDF.
  // For example, you can host it on Google Drive or your personal website.
  resume: '#'
};

export const PROJECTS: Project[] = [
  {
    id: 'project-one',
    title: 'Project One',
    description: 'A cutting-edge web application that revolutionizes task management for teams, built with React and a serverless backend.',
    detailedDescription: 'Project One was conceived to address the common pain points in team collaboration and task tracking. The core idea was to create a highly intuitive, real-time, and scalable platform where teams could seamlessly manage their workflows. We focused on a minimalist UI to reduce clutter and cognitive load, while powerful features like collaborative editing, smart notifications, and integration with popular developer tools run under the hood. The serverless architecture ensures high availability and cost-efficiency, scaling automatically with user demand.',
    imageUrl: 'https://picsum.photos/seed/project1/400/300',
    projectUrl: 'https://github.com/your-username-here/project-one',
    tags: ['React', 'TypeScript', 'Serverless', 'WebSockets']
  },
  {
    id: 'project-two',
    title: 'Project Two',
    description: 'An innovative e-commerce platform with a focus on user experience and performance, featuring a custom recommendation engine.',
    detailedDescription: 'The vision for Project Two was to build an e-commerce experience that feels personal and incredibly fast. Dissatisfied with generic, slow-loading online stores, we engineered a platform from the ground up using Next.js for server-side rendering and static site generation, achieving sub-second page loads. The centerpiece is a custom-built recommendation engine that uses machine learning to analyze user behavior and suggest products they\'ll genuinely love, moving beyond simple "customers who bought this also bought" logic. This creates a more engaging and effective shopping journey.',
    imageUrl: 'https://picsum.photos/seed/project2/400/300',
    projectUrl: 'https://github.com/your-username-here/project-two',
    tags: ['Next.js', 'GraphQL', 'Stripe', 'Machine Learning']
  },
  {
    id: 'project-three',
    title: 'Project Three',
    description: 'A data visualization dashboard that provides real-time insights from complex datasets, using D3.js for interactive charts.',
    detailedDescription: 'Project Three was born from the need to make complex data accessible and understandable to non-technical stakeholders. The idea was to transform raw, overwhelming datasets into beautiful, interactive stories. Using D3.js, we created a library of bespoke, animated charts and graphs that allow users to drill down into the data, uncover trends, and export insights. The dashboard connects to multiple data sources in real-time via WebSockets, ensuring that the visualizations always reflect the latest information, enabling data-driven decisions on the fly.',
    imageUrl: 'https://picsum.photos/seed/project3/400/300',
    projectUrl: 'https://github.com/your-username-here/project-three',
    tags: ['D3.js', 'React', 'WebSocket', 'Data Visualization']
  },
  {
    id: 'project-four',
    title: 'Project Four',
    description: 'A mobile-first social networking app designed to connect people with shared interests in their local area.',
    detailedDescription: 'In an increasingly digital world, Project Four aims to foster real-world connections. The core concept is a "hyper-local" social network that helps people discover and connect with others who share their hobbies and passions within their immediate vicinity. Built with React Native for a consistent cross-platform experience, the app features event creation, group chats, and a map-based discovery view. We prioritized privacy and safety, implementing robust user verification and content moderation systems to create a trustworthy community space.',
    imageUrl: 'https://picsum.photos/seed/project4/400/300',
    projectUrl: 'https://github.com/your-username-here/project-four',
    tags: ['React Native', 'Firebase', 'Redux', 'Geolocation']
  },
  {
    id: 'project-five',
    title: 'Project Five',
    description: 'A developer tool that automates code quality checks and deployment pipelines, improving team productivity.',
    detailedDescription: 'Project Five is a tool built by developers, for developers. The idea was to streamline the development lifecycle by automating tedious but critical tasks. It integrates directly with Git repositories to run a configurable pipeline of code linting, unit testing, integration testing, and vulnerability scanning on every commit. If all checks pass, it automatically handles the build and deployment process to various environments (staging, production). This CI/CD tool provides a centralized dashboard for viewing build statuses and logs, significantly reducing manual overhead and improving code quality.',
    imageUrl: 'https://picsum.photos/seed/project5/400/300',
    projectUrl: 'https://github.com/your-username-here/project-five',
    tags: ['Node.js', 'Docker', 'CI/CD', 'GitHub API']
  },
  {
    id: 'project-six',
    title: 'Project Six',
    description: 'An open-source component library for React, offering a set of accessible and themeable UI elements for developers.',
    detailedDescription: 'The goal of Project Six was to accelerate frontend development without sacrificing quality or accessibility. We developed an open-source library of common UI components (buttons, modals, forms, etc.) that are fully compliant with WCAG accessibility standards out-of-the-box. The library is built with a theming API, allowing developers to easily customize the look and feel to match their brand. Documented extensively with Storybook, it serves as a reliable foundation for building consistent and high-quality user interfaces, allowing teams to focus on application logic rather than reinventing UI primitives.',
    imageUrl: 'https://picsum.photos/seed/project6/400/300',
    projectUrl: 'https://github.com/your-username-here/project-six',
    tags: ['React', 'Storybook', 'Styled Components', 'Accessibility']
  }
];