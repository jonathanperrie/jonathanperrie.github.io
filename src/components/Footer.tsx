import React from 'react';

const Footer: React.FC = () => {
  return (
    <footer className="bg-green-800 border-t border-green-700 mt-20">
      <div className="container mx-auto px-4 py-8 flex justify-center items-center">
        <p className="text-green-400">
          &copy; {new Date().getFullYear()} Hello. All rights reserved.
        </p>
      </div>
    </footer>
  );
};

export default Footer;
