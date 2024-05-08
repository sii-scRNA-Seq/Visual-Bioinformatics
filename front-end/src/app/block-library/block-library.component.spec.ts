import { By } from '@angular/platform-browser';
import { ComponentFixture, TestBed } from '@angular/core/testing';
import { HttpClientTestingModule } from '@angular/common/http/testing';
import { MatCardModule } from '@angular/material/card';
import { MatSnackBarModule } from '@angular/material/snack-bar';

import { BlockLibraryComponent } from './block-library.component';
import { BlockService } from '../block.service';
import { MockBlockService } from '../mock-block.service';
import { MockOutputService } from '../mock-output.service';
import { OutputService } from '../output.service';

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
      providers: [
        { provide: BlockService, useClass: MockBlockService },
        { provide: OutputService, useClass: MockOutputService },
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

    it ('should call blockService.addBlock with pca when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const button = fixture.debugElement.query(By.css('#pca'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.addBlock).toHaveBeenCalledOnceWith('pca');
    });

    it ('should call blockService.addBlock with runumap when add button is clicked', () => {
      const blockService: BlockService = TestBed.inject(BlockService);
      spyOn(blockService, 'addBlock');
      const button = fixture.debugElement.query(By.css('#run-umap'));
      button.triggerEventHandler('click', {});
      fixture.detectChanges();
      expect(blockService.addBlock).toHaveBeenCalledOnceWith('runumap');
    });

    it ('should be available when blocks are not being executed', () => {
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('#load-data')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#basic-filtering')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#qc-plots')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#qc-filtering')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#variable-genes')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#pca')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#run-umap')).nativeElement.disabled).toEqual(false);
    });

    it ('should become disabled while blocks are being executed', () => {
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('#load-data')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#basic-filtering')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#qc-plots')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#qc-filtering')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#variable-genes')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#pca')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#run-umap')).nativeElement.disabled).toEqual(false);
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('#load-data')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#basic-filtering')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#qc-plots')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#qc-filtering')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#variable-genes')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#pca')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#run-umap')).nativeElement.disabled).toEqual(true);
    });

    it ('should become available once blocks have stopped being executed', () => {
      component.executingBlocks = true;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('#load-data')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#basic-filtering')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#qc-plots')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#qc-filtering')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#variable-genes')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#pca')).nativeElement.disabled).toEqual(true);
      expect(fixture.debugElement.query(By.css('#run-umap')).nativeElement.disabled).toEqual(true);
      component.executingBlocks = false;
      fixture.detectChanges(); 
      expect(fixture.debugElement.query(By.css('#load-data')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#basic-filtering')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#qc-plots')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#qc-filtering')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#variable-genes')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#pca')).nativeElement.disabled).toEqual(false);
      expect(fixture.debugElement.query(By.css('#run-umap')).nativeElement.disabled).toEqual(false);
    });
  });
});
