import React from 'react';
import { Project } from '../types';

interface ProjectCardProps {
  project: Project;
  onNavigate: (page: 'project', projectId: string) => void;
}

const ProjectCard: React.FC<ProjectCardProps> = ({ project, onNavigate }) => {
  return (
    <div 
      className="bg-green-950 rounded-lg overflow-hidden shadow-lg hover:shadow-lime-500/50 transition-shadow duration-300 flex flex-col cursor-pointer"
      onClick={() => onNavigate('project', project.id)}
    >
      <img src={project.imageUrl} alt={project.title} className="w-full h-48 object-cover" />
      <div className="p-6 flex flex-col flex-grow">
        <h3 className="text-xl font-bold text-green-200 mb-2">{project.title}</h3>
        <p className="text-green-400 mb-4 flex-grow">{project.description}</p>
        <div className="mb-4">
          {project.tags.map((tag, index) => (
            <span key={index} className="inline-block bg-lime-900/50 text-lime-300 text-xs font-semibold mr-2 mb-2 px-2.5 py-0.5 rounded-full">
              {tag}
            </span>
          ))}
        </div>
        <div className="mt-auto">
          <span
            className="text-lime-400 group-hover:text-lime-300 font-semibold transition-colors"
          >
            Learn More
          </span>
        </div>
      </div>
    </div>
  );
};

export default ProjectCard;
