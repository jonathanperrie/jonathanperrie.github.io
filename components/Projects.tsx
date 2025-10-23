
import React from 'react';
import { PROJECTS } from '../constants';
import ProjectCard from './ProjectCard';

interface ProjectsProps {
  onNavigate: (page: 'project', projectId: string) => void;
}

const Projects: React.FC<ProjectsProps> = ({ onNavigate }) => {
  return (
    <section id="projects" className="py-20 md:py-32">
      <h2 className="text-3xl md:text-4xl font-bold text-center text-green-100 mb-16">My Projects</h2>
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
        {PROJECTS.map((project, index) => (
          <ProjectCard key={index} project={project} onNavigate={onNavigate} />
        ))}
      </div>
    </section>
  );
};

export default Projects;