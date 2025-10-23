
import React from 'react';

interface IconProps {
  className?: string;
}

const GoogleScholarIcon: React.FC<IconProps> = ({ className }) => (
  <svg
    xmlns="http://www.w3.org/2000/svg"
    viewBox="0 0 24 24"
    fill="currentColor"
    className={className}
  >
    <path d="M12 24a7 7 0 1 1 0-14 7 7 0 0 1 0 14zm0-24L0 6v5.217l12 6.817 12-6.817V6L12 0z" />
  </svg>
);

export default GoogleScholarIcon;
