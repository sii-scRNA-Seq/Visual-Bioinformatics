import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { BlockLibraryComponent } from './block-library.component';
import { BlockService } from '../block.service';

describe('BlockLibraryComponent', () => {
  let component: BlockLibraryComponent;
  let fixture: ComponentFixture<BlockLibraryComponent>;

  beforeEach(() => {
    TestBed.configureTestingModule({
      declarations: [BlockLibraryComponent],
      imports: [
        HttpClientTestingModule,
        MatCardModule,
        MatSnackBarModule,
      ],
    });
    fixture = TestBed.createComponent(BlockLibraryComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  describe('Adding Blocks From Library', () => {
    it ('should call blockService.addBlock with loaddata when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const button = fixture.debugElement.query(By.css('#load-data'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.addBlock).toHaveBeenCalledOnceWith('loaddata');
    });

    it ('should call blockService.addBlock with basicfiltering when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const button = fixture.debugElement.query(By.css('#basic-filtering'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.addBlock).toHaveBeenCalledOnceWith('basicfiltering');
    });

    it ('should call blockService.addBlock with qcplots when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const button = fixture.debugElement.query(By.css('#qc-plots'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.addBlock).toHaveBeenCalledOnceWith('qcplots');
    });

    it ('should call blockService.addBlock with qcfiltering when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const button = fixture.debugElement.query(By.css('#qc-filtering'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.addBlock).toHaveBeenCalledOnceWith('qcfiltering');
    });

    it ('should call blockService.addBlock with variablegenes when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const button = fixture.debugElement.query(By.css('#variable-genes'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.addBlock).toHaveBeenCalledOnceWith('variablegenes');
    });
  });
});
